#!/usr/bin/env node

"use strict";

const EventEmitter = require('events');
const assert       = require('assert');
const LOGGER       = require('../lib/logsetup.js');

const config       = require('../lib/config.js');
const RangerServer = require('../lib/server.js');

/**
 * NodeJS implementation of cluster
 * @external cluster
 * @see  {@link https://nodejs.org/dist/latest-v4.x/docs/api/cluster.html|cluster}
 */

/**
 * <p>Executable server module</p>
 *
 * <p>If this module is called as executable it will run a Ranger Server with
 * available configuration.</p>
 *
 * <p>Providing config settings</br>
 * Settings are provided from 3 locations:<p/>
 * <ol>
 *   <li>Command line - run with -h to see options.</li>
 *   <li>Config json file can be read if it is provided on command line</br>
 *       by running with -c PATH or --configfile=PATH</li>
 *   <li>There are some defaults, which can be found in lib/config.js</li>
 * </ol>
 *
 * Call to <em>config.provide()</em> must occur here before requiring
 * controller so that options object is built before it is provided to the
 * other modules
 *
 * @module bin/server
 * @requires {@link external:events|events}
 * @requires {@link external:assert|assert}
 * @requires {@link external:cluster|cluster}
 *
 * @author Andrew Nowak
 * @author Jaime Tovar
 * @copyright Genome Research Limited 2016
 */

/**
* Event emitted when a server starts. Usually bubbling up from RangerServer.
* @type {String}
*/
const SERVER_STARTED        = RangerServer.SERVER_STARTED;
/**
* Event emitted when a server closes. Usually bubbling up from RangerServer.
* @type {String}
*/
const SERVER_CLOSED         = RangerServer.SERVER_CLOSED;
/**
* Emitted when cluster starts
* @type {String}
*/
const CLUSTER_STARTED       = 'cluster_started';
/**
* Emitted each time the worker starts
* @type {String}
*/
const WORKER_STARTED        = 'worker_started';
/**
* Emitted each time the cluster forks a new worker
* @type {String}
*/
const WORKER_FORKED         = 'worker_forked';
/**
* Emitted when the cluster notices a worker has gone down
* @type {String}
*/
const WORKER_CLOSED         = 'worker_closed';
/**
* Emitted when the cluster identifies too many workers has died in a short span
* of time. The cluster will stop forking after this event is emitted. Short
* after it will try to clean up and exit.
* @type {String}
*/
const HIT_LIMIT_CONSEC_FORK = 'hit_limit_consec_forks';
/**
* Set to <i>210</i>. Error code at process exit when too many forks died in
* short span of time.
* @type {Number}
*/
const ERROR_SERVER_LIMIT_CONSEC_FORK = 210;


/**
 * Parent class for brokers
 */
class Broker extends EventEmitter {
  constructor(serverFactory) {
    super();
    this.serverFactory = serverFactory;
    this.serverFactory.on(SERVER_STARTED, ( server ) => {
      this.emit(SERVER_STARTED, server);
    });
    this.serverFactory.on(SERVER_CLOSED, ( server ) => {
      this.emit(SERVER_CLOSED, server);
    });
  }

  start() {
    LOGGER.debug('broker start server');
  }
}

/**
 * Simplest broker, will immediatelly delegate to the server factory to start
 * a server.
 */
class FlatBroker extends Broker {
  start() {
    super.start();
    this.serverFactory.startServer();
  }
}

/**
 * Broker which uses the {@link external:cluster|cluster module} to orchestrate
 * a cluster of servers.
 */
class ClusteredBroker extends Broker {
  /**
   * Creates a cluster and runs multiple workers. Emits events when the cluster
   * starts and when workers are forked, when they start and die. If too many wokers
   * die in short time span, the master will stop forking workers, clean-up and exit.
   */
  start() {
    super.start();
    let options = config.provide();
    let numworkers = options.get('numworkers');
    assert(Number.isInteger(numworkers), 'numworkers must be an integer');
    const cluster = require('cluster');
    let self = this;
    if ( cluster.isMaster ) {
      LOGGER.info(config.logOpts());
      for (let i = 0; i < numworkers; i++) {
        cluster.fork();
        self.emit(WORKER_FORKED);
      }
      let consec = 0;
      let exiting = false;
      cluster.on('exit', (worker, code, signal) => {
        LOGGER.debug('Worker %d died (%s). Forking to replace ...', worker.id, signal || code);
        let waitingConsec = options.get('clustertimeout');
        let maxConsec = options.get('clustermaxdeaths');
        LOGGER.debug(`${consec} forks have died in the previous ${waitingConsec} seconds.`);
        if ( consec >= maxConsec ) {
          if ( !exiting ) {
            exiting = true;
            LOGGER.error('Too many forks started in short span of time. Trying to exit now.');
            self.emit(HIT_LIMIT_CONSEC_FORK);
            setTimeout( () => {
              process.exit(ERROR_SERVER_LIMIT_CONSEC_FORK);
            }, 3000);
          }
        } else {
          consec += 1;
          cluster.fork();
          self.emit(WORKER_FORKED);
          setTimeout( () => {
            consec -= 1;
          }, waitingConsec * 1000 );
        }
      });
      self.emit(CLUSTER_STARTED, cluster);
    } else {
      if ( cluster.isWorker ) {
        self.emit(WORKER_STARTED, cluster.worker);
        LOGGER.debug('WORKER: new fork ' + cluster.worker.id);
        self.serverFactory.startServer();
      }
    }
  }
}

/**
 * Factory for brokers
 *
 * @example
 *
 * const RangerServer = require('../lib/server.js');
 *
 * let brokerFactory = new BrokerFactory();
 * let serverFactory = new RangerServer.ServerFactory();
 * let broker = brokerFactory.buildBroker(serverFactory);
 *
 * process.nextTick(() => {
 *   broker.start();
 * });
 */
class BrokerFactory {
  /**
   * Use the number of workers build the specific broker needed. Pass the server
   * factory to the broker so it can build servers with it.
   */
  buildBroker(serverFactory) {
    assert(serverFactory, 'serverFactory is required');
    let options = config.provide();
    let numworkers = options.get('numworkers');
    assert(Number.isInteger(numworkers), 'numworkers must be an integer');
    let broker;
    if ( !numworkers ) {
      broker = new FlatBroker(serverFactory);
    } else {
      broker = new ClusteredBroker(serverFactory);
    }
    return broker;
  }
}

/**
 * Application's main method.
 */
if ( require.main === module ) {
  const options = config.provide(config.fromCommandLine);

  if ( options.get('debug') ) {
    LOGGER.level = 'debug';
  }

  let bf = new BrokerFactory();
  let sf = new RangerServer.ServerFactory();

  sf.on(SERVER_STARTED, () => { LOGGER.debug('Server factory started server'); });
  sf.on(SERVER_CLOSED,  () => { LOGGER.debug('Server factory server closed'); });

  let broker = bf.buildBroker(sf);

  broker.on(CLUSTER_STARTED, () => { LOGGER.debug('Cluster started'); });
  broker.on(WORKER_STARTED,  () => { LOGGER.debug('Cluster - worker started'); });
  broker.on(WORKER_FORKED,   () => { LOGGER.debug('Cluster - worker forked'); });

  process.nextTick(() => {
    broker.start();
  });
}

module.exports = {
  BrokerFactory:   BrokerFactory,
  // Event names
  CLUSTER_STARTED:       CLUSTER_STARTED,
  WORKER_STARTED:        WORKER_STARTED,
  WORKER_FORKED:         WORKER_FORKED,
  WORKER_CLOSED:         WORKER_CLOSED,
  SERVER_STARTED:        SERVER_STARTED,
  SERVER_CLOSED:         SERVER_CLOSED,
  HIT_LIMIT_CONSEC_FORK: HIT_LIMIT_CONSEC_FORK
};
