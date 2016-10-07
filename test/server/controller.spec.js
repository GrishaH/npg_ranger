/* globals describe, it, expect, beforeAll, afterAll */

"use strict";

const assert  = require('assert');
const http    = require('http');
const fs      = require('fs');
const fse     = require('fs-extra');
const tmp     = require('tmp');
const RangerController = require('../../lib/server/controller.js');
const config  = require('../../lib/config.js');

// Create temp dir here so it is available for all tests.
// Use this dir as a default dir that will be available in all.
var tmpDir    = config.tempFilePath('npg_ranger_controller_test_');
var dummy     = function() { return {tempdir: tmpDir, skipauth: true}; };
var options;

describe('Creating object instance - synch', function() {
  beforeAll(function() {
    options = config.provide(dummy);
    fse.ensureDirSync(tmpDir);
  });

  afterAll(function() {
    fse.removeSync(tmpDir);
  });

  it('request object is not given - error', function() {
    expect( () => {new RangerController();} ).toThrowError(assert.AssertionError,
    'HTTP request object is required');
  });
  it('request is not an object - error', function() {
    expect( () => {new RangerController(1);} ).toThrowError(assert.AssertionError,
    'HTTP request object is required');
  });
  it('response object is not given - error', function() {
    expect( () => {new RangerController({});} ).toThrowError(assert.AssertionError,
    'Server response object is required');
  });
  it('response is not an object - error', function() {
    expect( () => {new RangerController({}, 1);} ).toThrowError(assert.AssertionError,
    'Server response object is required');
  });
});

describe('set error response', function() {
  // Create server HTTP server object.
  const server = http.createServer();
  // Generate synchronously a temporary file name.
  var socket = tmp.tmpNameSync();

  beforeAll((done) => {
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
    // Start listening on a socket
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });

  // This tidy-up callback is not called when the spec exits early
  // due to an error. Seems to be a bug in jasmine.
  afterAll(function() {
    server.close();
    try { fs.unlinkSync(socket); } catch ( e ) { console.log(e); }
    fse.removeSync(tmpDir);
  });

  it('db object is not given or is not an object - error', function(done) {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      assert(typeof request == 'object');
      expect( () => {new RangerController(request, response);} ).toThrowError(
        assert.AssertionError, 'DB handle object is required');
      expect( () => {new RangerController(request, response, 1);} ).toThrowError(
        assert.AssertionError, 'DB handle object is required');
      response.write('payload');
      response.end();
    });

    http.get({socketPath: socket}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        done();
      });
    });
  });

  it('Setting values of instance variables', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c;
      expect( () => {c = new RangerController(request, response, {one: "two"});} ).not.toThrow();
      expect((typeof c == 'object')).toBe(true);
      expect((c instanceof RangerController)).toBe(true);
      expect((c.request == request)).toBe(true);
      expect((c.response == response)).toBe(true);
      expect(c.db).toEqual({one: "two"});
      response.end();
      done();
    });
    http.get({socketPath: socket}, function() {});
  });
});

describe('Handling requests - error responses', function() {
  const server = http.createServer();
  var socket = tmp.tmpNameSync();
  beforeAll((done) => {
    options = config.provide(dummy);
    fse.ensureDirSync(tmpDir);
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });

  afterAll(function() {
    server.close();
    try { fs.unlinkSync(socket); } catch (e) {}
    fse.removeSync(tmpDir);
  });

  it('Method not allowed error', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });

    let options = {socketPath: socket, method:    'POST'};
    let req = http.request(options);
    req.on('response', (response) => {
      expect(response.headers['content-type']).toEqual('application/json');
      expect(response.statusCode).toEqual(405);
      expect(response.statusMessage).toEqual('POST request is not allowed');
      done();
    });
    req.end();
  });

  it('Authentication error', function(done) {
    config.provide( () => {
      return {tempdir:  tmpDir,
              skipauth: false};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });
    http.get({socketPath: socket}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(401);
        expect(response.statusMessage).toEqual('Proxy authentication required');
        expect(JSON.parse(body)).toEqual(
          {error: {type:    "InvalidAuthentication",
                   message: "Proxy authentication required"}});
        done();
      });
    });
  });

  it('Not found error, no auth', function(done) {
    config.provide( () => {
      return {tempdir:  tmpDir,
              skipauth: true};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });
    http.get({socketPath: socket, path: '/invalid'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual('URL not found : /invalid');
        expect(JSON.parse(body)).toEqual(
          {error: {type:    "NotFound",
                   message: "URL not found : /invalid"}});
        done();
      });
    });
  });

  it('Not found error, auth checked', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });
    let req = http.request({socketPath: socket, path: '/invalid'});
    req.setHeader('X-Remote-User', 'user1');
    req.on('response', function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual('URL not found : /invalid');
        expect(JSON.parse(body)).toEqual(
          {error: {type:    "NotFound",
                   message: "URL not found : /invalid"}});
        done();
      });
    });
    req.end();
  });


  it('Invalid input error for a sample url', ( done ) => {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });

    http.get({socketPath: socket, path: '/sample'}, ( response ) => {
      var body = '';
      response.on('data', ( d ) => { body += d;});
      response.on('end', () => {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(422);
        let m = 'Invalid request: sample accession number should be given';
        expect(response.statusMessage).toEqual(m);
        expect(JSON.parse(body)).toEqual(
          {error: {type:    "InvalidInput",
                   message: m}});
        done();
      });
    });
  });

  it('Invalid input error for a file url', ( done ) => {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
    });

    http.get({socketPath: socket, path: '/file'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(422);
        let m = 'Invalid request: file name should be given';
        expect(response.statusMessage).toEqual(m);
        expect(JSON.parse(body)).toEqual(
          {error: {type:    "InvalidInput",
                   message: m}});
        done();
      });
    });
  });

  it('Invalid input error for a vcf file when multiref set', (done) => {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      // Set multiref mode
      config.provide(() => { return {tempdir: tmpDir, multiref: true, skipauth: true}; });
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).not.toThrow();
      // unset multiref
      config.provide(dummy);
    });

    http.get({socketPath: socket, path: '/sample?accession=XYZ120923&format=vcf'}, (response) => {
      var body = '';
      response.on('data', (d) => { body += d;});
      response.on('end', () => {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(422);
        let m = 'Invalid request: cannot produce VCF files while multiref set on server';
        expect(response.statusMessage).toEqual(m);
        expect(JSON.parse(body)).toEqual(
          {
            error: {
              type: "InvalidInput",
              message: m
            }
          });
        done();
      });
    });
  });
});

describe('Redirection in json response', function() {
  const server = http.createServer();
  var socket = tmp.tmpNameSync();
  let id          = 'EGA45678';
  let server_path_basic = '/ga4gh/v.0.1/get/sample';
  let server_path = server_path_basic + '/' + id;

  beforeAll((done) =>  {
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {});
      c.handleRequest();
    });
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });

  afterAll(function() {
    server.close();
    try { fs.unlinkSync(socket); } catch (e) {}
    fse.removeSync(tmpDir);
  });

  it('invalid url - no id - error response', function(done) {
    http.get({socketPath: socket, path: server_path_basic}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual('URL not found : ' + server_path_basic);
        done();
      });
    });
  });

  it('invalid url - no id - error response', function(done) {
    let path = server_path_basic + '/';
    http.get({socketPath: socket, path: path}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual('URL not found : ' + path);
        done();
      });
    });
  });

  it('invalid sample id - error response', function(done) {
    let path = server_path_basic + 'ERS-4556';
    http.get({socketPath: socket, path: path}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(404);
        expect(response.statusMessage).toEqual('URL not found : ' + path);
        done();
      });
    });
  });

  it('successful redirection, no query params', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM`;
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  it('successful redirection, format given', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?format=CRAM'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=CRAM`;
        expect(JSON.parse(body)).toEqual({format: 'CRAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  it('successful redirection, chromosome given', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1`;
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  it('successful redirection, range start given', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&start=3'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1%3A4`;
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  it('successful redirection, range end given', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&end=4'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1%3A1-4`;
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  it('successful redirection, range start and end given', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&start=4&end=400'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(200);
        expect(response.statusMessage).toEqual(
          'OK, see redirection instructions in the body of the message');
        let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1%3A5-400`;
        expect(JSON.parse(body)).toEqual({format: 'BAM', urls: [{'url': url}]});
        done();
      });
    });
  });

  ['bam', 'BAM', 'sam', 'SAM', 'cram', 'CRAM', 'vcf', 'VCF'].forEach( ( value ) => {
    it('successful redirection, query with all possible params', function(done) {
      http.get(
        { socketPath: socket,
          path: server_path + `?referenceName=chr1&start=4&end=400&format=${value}`}, function(response) {
        var body = '';
        response.on('data', function(d) { body += d;});
        response.on('end', function() {
          expect(response.headers['content-type']).toEqual('application/json');
          expect(response.statusCode).toBe(200);
          expect(response.statusMessage).toBe(
            'OK, see redirection instructions in the body of the message');
          let formatUpperCase = value.toUpperCase();
          let url = `http://localhost/sample?accession=${id}&format=${formatUpperCase}&region=chr1%3A5-400`;
          expect(JSON.parse(body)).toEqual({format: `${formatUpperCase}`, urls: [{'url': url}]});
          done();
        });
      });
    });
  });

  it('redirection error, range is given, reference is missing', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?start=4&end=400'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(422);
        expect(response.statusMessage).toBe(
          "'referenceName' attribute requered if 'start' or 'end' attribute is given");
        done();
      });
    });
  });

  it('redirection error, range start is not an integer', function(done) {
    http.get(
     {socketPath: socket,
      path: server_path + '?referenceName=chr1&start=5.5&end=400'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(422);
        expect(response.statusMessage).toEqual(
          "'5.5' is not an integer");
        done();
      });
    });
  });

  it('redirection error, range start is a negative integer', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&start=-44&end=400'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(422);
        expect(response.statusMessage).toEqual("'-44' is not an unsigned integer");
        done();
      });
    });
  });

  it('redirection error, range end is not an integer', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&start=4&end=foo'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(422);
        expect(response.statusMessage).toEqual("'foo' is not an integer");
        done();
      });
    });
  });

  it('redirection error, range end is a negative integer', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&start=4&end=-400'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(422);
        expect(response.statusMessage).toEqual("'-400' is not an unsigned integer");
        done();
      });
    });
  });

  it('redirection error, range start is bigger than range end', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?referenceName=chr1&start=400&end=4'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(422);
        expect(response.statusMessage).toEqual(
          'Range end should be bigger than start');
        done();
      });
    });
  });

  it('redirection error, unknown format requested', function(done) {
    http.get(
      { socketPath: socket,
        path: server_path + '?format=fa'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(409);
        expect(response.statusMessage).toEqual(
          "Format 'fa' is not supported, supported formats: BAM, CRAM, SAM, VCF");
        done();
      });
    });
  });

});

describe('content type', function() {
  const server = http.createServer();
  var socket = tmp.tmpNameSync();

  beforeAll((done) => {
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });
  afterAll(function() {
    server.close();
    try { fs.unlinkSync(socket); } catch (e) {}
    fse.removeSync(tmpDir);
  });

  it('data format driven content type', function(done) {
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.contentType();} )
        .toThrowError(assert.AssertionError,
        'Non-empty format string should be given');
      expect(c.contentType('SAM')).toBe('text/vnd.ga4gh.sam');
      expect(c.contentType('VCF')).toBe('text/vnd.ga4gh.vcf');
      expect(c.contentType('BAM')).toBe('application/vnd.ga4gh.bam');
      expect(c.contentType('CRAM')).toBe('application/vnd.ga4gh.cram');
      done();
    });

    http.get({socketPath: socket, path: '/file'}, function(response) {
      let body = '';
      response.on('data', function(d) { body += d;});
    });
  });
});

describe('trailers in response', function() {
  const server = http.createServer();
  var socket = tmp.tmpNameSync();

  beforeAll((done) => {
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
      done();
    });
  });
  afterAll(function() {
    server.close();
    try { fs.unlinkSync(socket); } catch (e) {}
    fse.removeSync(tmpDir);
  });

  it('no trailers without TE header', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect(c.sendTrailer).toBe(false);
      done();
    });

    http.get({socketPath:socket, path: '/file'}, function() {});
  });

  ['TE', 'te', 'Te', 'tE'].forEach( ( headerName ) => {
    it(`trailers with ${headerName} header`, function(done) {
      server.removeAllListeners('request');
      server.on('request', (request, response) => {
        let c = new RangerController(request, response, {one: "two"});
        expect(c.sendTrailer).toBe(true);
        done();
      });

      let headers = {};
      headers[headerName] = 'trailers';
      http.get({socketPath:socket, path: '/file', headers: headers}, () => {});
    });
  });
});

describe('CORS in response', function() {
  var server;
  let serverPath = '/ga4gh/v.0.1/get/sample/EGA45678';
  let socket = tmp.tmpNameSync();

  let checkHeaders = (headers, origin) => {
    expect(headers.vary).toBe('Origin', 'Vary header is set');
    expect(headers['access-control-allow-origin']).toBe(
      origin, `allowed origin is ${origin}`);
    expect(headers['access-control-allow-methods']).toBe(
      'GET,OPTIONS', 'allowed methods are set');
    expect(headers['access-control-allow-headers']).toBe(
      'TE,X-Remote-User', 'allowed headers are set');
    expect(headers['access-control-max-age']).toBe('1800', 'max age is set');
    expect(Object.keys(headers).indexOf('Access-Control-Allow-Credentials')).toBe(
      -1, 'Access-Control-Allow-Credentials header is not set');
  };

  beforeAll((done) => {
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
    server = http.createServer();
    server.listen(socket, () => { console.log('listening'); done();});
  });
  afterAll(function() {
    server.close();
    try { fs.unlinkSync(socket); } catch (e) {}
    fse.removeSync(tmpDir);
  });

  it('no CORS in a response to a standart request', function(done) {
    config.provide( () => {
      return {tempdir: tmpDir, anyorigin: false, originlist: null, skipauth: true};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect('origin' in request.headers).toBe(false, 'request does not have Origin header');
      let c = new RangerController(request, response, {one: "two"});
      c.handleRequest();
    });

    http.get({socketPath:socket, path: serverPath}, (res) => {
      expect( Object.keys(res.headers).filter((headerName) => {
        return headerName.startsWith('Access-Control');
      }).length).toBe(0, 'no CORS headers in reply');
      expect(res.headers.vary).toBe('Origin', 'Vary header is set');
      done();
    });
  });

  it('no CORS headers in a response to CORS GET request due to server options', function(done) {
    config.provide( () => {
      return {tempdir: tmpDir, anyorigin: false, originlist: null, skipauth: true};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect('origin' in request.headers).toBe(true, 'request has Origin header');
      let c = new RangerController(request, response, {one: "two"});
      c.handleRequest();
    });

    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://some.com'},
                   method:     'GET'};
    let req = http.request(options);
    req.on('response', (res) => {
      expect( Object.keys(res.headers).filter((headerName) => {
        return headerName.startsWith('access-control');
      }).length).toBe(0, 'no CORS headers in reply');
      expect(res.headers.vary).toBe('Origin', 'Vary header is set');
      done();
    });
    req.end();
  });

  it('no CORS headers in a response to CORS OPTIONS request due to server options', function(done) {
    config.provide( () => {
      return {tempdir: tmpDir, anyorigin: false, originlist: null, skipauth: true};
    });
    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://some.com'},
                   method:     'OPTIONS'};
    let req = http.request(options);
    req.on('response', (res) => {
      expect( Object.keys(res.headers).filter((headerName) => {
        return headerName.startsWith('access-control');
      }).length).toBe(0, 'no CORS headers in reply');
      done();
    });
    req.end();
  });

  it('Allow all CORS headers in a response to CORS GET request', function(done) {
    config.provide( () => {
      return {tempdir: tmpDir, anyorigin: true, originlist: null, skipauth: true};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect('origin' in request.headers).toBe(true, 'request has Origin header');
      let c = new RangerController(request, response, {one: "two"});
      c.handleRequest();
    });

    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://some.com'},
                   method:     'GET'};
    let req = http.request(options);
    req.on('response', (res) => {
      checkHeaders(res.headers, '*');
      done();
    });
    req.end();
  });

  it('Allow all CORS headers in a response to CORS OPTIONS request', function(done) {
    config.provide( () => {
      return {tempdir: tmpDir, anyorigin: true, originlist: null, skipauth: true};
    });
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect('origin' in request.headers).toBe(true, 'request has Origin header');
      let c = new RangerController(request, response, {one: "two"});
      c.handleRequest();
    });

    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://some.com'},
                   method:     'OPTIONS'};
    let req = http.request(options);
    req.on('response', (res) => {
      checkHeaders(res.headers, '*');
      done();
    });
    req.end();
  });

  it('No CORS headers since the origin is not white listed', function(done) {
    config.provide( () => {
      return {tempdir:     tmpDir,
              anyorigin:   false,
              originlist: 'http://other.com,http://other.com:9090',
              skipauth:   true};
    });
    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://some.com'},
                   method:     'OPTIONS'};
    let req = http.request(options);
    req.end();
    req.on('response', (res) => {
      expect( Object.keys(res.headers).filter((headerName) => {
        return headerName.startsWith('access-control');
      }).length).toBe(0, 'no CORS headers in reply');
      done();
    });
  });

  it('Origin-specific CORS headers set for a white listed origin', function(done) {
    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://other.com'},
                   method:     'OPTIONS'};
    let req = http.request(options);
    req.on('response', (res) => {
      checkHeaders(res.headers, 'http://other.com');
      done();
    });
    req.end();
  });

  it('Additional CORS header is set when running with authorization', function(done) {
    config.provide( () => {
      return {tempdir:    tmpDir,
              anyorigin:  false,
              originlist: 'http://other.com,http://other.com:9090',
              skipauth:   false};
    });
    let options = {socketPath: socket,
                   path:       serverPath,
                   headers:    {Origin: 'http://other.com:9090'},
                   method:    'OPTIONS'};
    let req = http.request(options);
    req.on('response', (res) => {
      expect(res.headers['access-control-allow-origin']).toBe(
        'http://other.com:9090', 'origin http://other.com:9090 is allowed');
      expect(res.headers['access-control-allow-credentials']).toBe(
        'true', 'Access-Control-Allow-Credentials header is set');
      done();
    });
    req.end();
  });
});