/* globals describe, it, expect, beforeAll, afterAll*/

"use strict";

const assert  = require('assert');
const http    = require('http');
const fs      = require('fs');
const fse     = require('fs-extra');
const os      = require('os');
const path    = require('path');
const tmp     = require('tmp');
const RangerController = require('../../lib/server/controller.js');
const config  = require('../../lib/config.js');

// Create temp dir here so it is available for all tests.
// Use this dir as a default dir that will be available in all.
var tmpDir    = config.tempFilePath('npg_ranger_controller_test_');
var dummy     = function() { return {tempdir: tmpDir}; };
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
  it('default configs are generated by controller constructor', function() {
    let c;
    expect( () => {c = new RangerController({}, {}, {});} ).not.toThrow();
    expect(c.tmpDir === options.get('tempdir')).toBe(true);
    expect(fs.existsSync(c.tmpDir));
    expect(c.skipAuth).toBe(false);
  });
  it('passing configs to provide overwrites defaults', function() {
    let temporaryDir = tmp.dirSync({ template: path.join(os.tmpdir(), 'tmp-XXXXXX') });
    options = config.provide( () => {
      return {tempdir: temporaryDir.name, skipauth: true};
    });
    let c;
    expect( () => {c = new RangerController({}, {}, {});} ).not.toThrow();
    expect(c.tmpDir === temporaryDir.name).toBe(true);
    expect(c.skipAuth).toBe(true);
    temporaryDir.removeCallback();
  });
  it('configs provided in constructor overwrite defaults', function() {
    let c;
    let temporaryDir = tmp.dirSync({ template: path.join(os.tmpdir(), 'tmp-XXXXXX') });
    expect( () => {c = new RangerController({}, {}, {}, temporaryDir.name, true);} ).not.toThrow();
    expect(c.tmpDir === temporaryDir.name).toBe(true);
    expect(c.skipAuth).toBe(true);
    expect( () => {c = new RangerController({}, {}, {}, temporaryDir.name, false);} ).not.toThrow();
    expect(c.tmpDir === temporaryDir.name).toBe(true);
    expect(c.skipAuth).toBe(false);
    temporaryDir.removeCallback();
  });
});

describe('set error response', function() {
  // Create server HTTP server object.
  const server = http.createServer();
  // Generate synchronously a temporary file name.
  var socket = tmp.tmpNameSync();

  beforeAll(function() {
    // Start listening on a socket
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
    });
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
  });

  // This tidy-up callback is not called when the spec exits early
  // due to an error. Seems to be a bug in jasmine.
  afterAll(function() {
    server.close();
    try { fs.unlinkSync(socket); } catch (e) {}
    fse.removeSync(tmpDir);
  });

  it('db object is not given or is not an object - error', function(done) {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      console.log(response.toString());
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
      expect(c.tmpDir).toBe(tmpDir);
      expect(c.skipAuth).toBe(false);
      expect(c.unsafe).toBe(false);
      expect( () => {c = new RangerController(request, response, {}, null, 0);} ).not.toThrow();
      expect(c.tmpDir).toBe(tmpDir);
      expect(c.skipAuth).toBe(false);
      expect(c.unsafe).toBe(false);
      expect( () => {c = new RangerController(request, response, {}, '', true);} ).not.toThrow();
      expect(c.tmpDir).toBe(tmpDir);
      expect(c.skipAuth).toBe(true);
      expect(c.unsafe).toBe(false);

      // Set unsafe mode
      config.provide(function() { return {tempdir: tmpDir, unsafe: true}; });
      expect( () => {c = new RangerController(request, response, {}, '', false);} ).not.toThrow();
      expect(c.skipAuth).toBe(false);
      expect(c.unsafe).toBe(true);
      expect( () => {c = new RangerController(request, response, {}, '', true);} ).not.toThrow();
      expect(c.skipAuth).toBe(true);
      expect(c.unsafe).toBe(true);

      // Reset safe mode
      config.provide(dummy);

      response.end();
      done();
    });
    http.get({socketPath: socket}, function() {});
  });

  it('Temporary directory should exist', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      expect( () => {new RangerController(request, response, {}, '/some/dir');} )
        .toThrowError(assert.AssertionError,
        "Temp data directory '/some/dir' does not exist");
      let c;
      expect( () => {c = new RangerController(request, response, {}, 'test', 0);} ).not.toThrow();
      expect(c.tmpDir).toBe('test');
      done();
    });
    http.get({socketPath: socket}, function() {});
  });
});

describe('Handling requests - error responses', function() {
  const server = http.createServer();
  var socket = tmp.tmpNameSync();
  beforeAll(function() {
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
    });
    options = config.provide(dummy);
    fse.ensureDirSync(tmpDir);
  });

  afterAll(function() {
    server.close();
    try { fs.unlinkSync(socket); } catch (e) {}
    fse.removeSync(tmpDir);
  });

  it('Data host name argument is required', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect( () => {c.handleRequest();} ).toThrowError(assert.AssertionError,
        'The data host name is required');
      c.handleRequest('localhost');
      response.end();
      done();
    });
    http.get({socketPath: socket}, function() {});
  });

  it('Authentication error', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"});
      expect(c.skipAuth).toBe(false);
      expect( () => {c.handleRequest('localhost');} ).not.toThrow();
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
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"}, null, true);
      expect(c.skipAuth).toBe(true);
      expect( () => {c.handleRequest('localhost');} ).not.toThrow();
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
      let c = new RangerController(request, response, {one: "two"}, null, false);
      expect(c.skipAuth).toBe(false);
      expect( () => {c.handleRequest('localhost');} ).not.toThrow();
    });
    let req = http.request({socketPath: socket, path: '/invalid'});
    req.setHeader('X-Remote-User', 'user1');
    req.end();
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
  });


  it('Invalid input error for a sample url', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"}, null, true);
      expect(c.skipAuth).toBe(true);
      expect( () => {c.handleRequest('localhost');} ).not.toThrow();
    });

    http.get({socketPath: socket, path: '/sample'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
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

  it('Invalid input error for a file url', function(done) {
    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"}, null, true);
      expect(c.skipAuth).toBe(true);
      expect( () => {c.handleRequest('localhost');} ).not.toThrow();
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

  it('Invalid input error for a vcf file when safe mode disabled', function(done) {

    server.removeAllListeners('request');
    server.on('request', (request, response) => {
      // Set unsafe mode
      config.provide(function() { return {tempdir: tmpDir, unsafe: true}; });

      let c = new RangerController(request, response, {one: "two"}, null, true);
      expect(c.skipAuth).toBe(true);
      expect(c.unsafe).toBe(true);
      expect( () => {c.handleRequest('localhost');} ).not.toThrow();

      // Reset safe mode
      config.provide(dummy);
    });

    http.get({socketPath: socket, path: '/sample?accession=XYZ120923&format=vcf'}, function(response) {
      var body = '';
      response.on('data', function(d) { body += d;});
      response.on('end', function() {
        expect(response.headers['content-type']).toEqual('application/json');
        expect(response.statusCode).toEqual(422);
        let m = 'Invalid request: cannot produce VCF files while server is not in safe mode';
        expect(response.statusMessage).toEqual(m);
        expect(JSON.parse(body)).toEqual(
          {error: {type:    "InvalidInput",
                   message: m}});
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

  beforeAll(function() {
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
    });
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {}, null, true);
      c.handleRequest('localhost');
    });
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
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
        //path: server_path}, function(response) {
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
        let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1%3A1-5`;
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
        let url = `http://localhost/sample?accession=${id}&format=BAM&region=chr1%3A5-401`;
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
          let formatUpperCase = value.toUpperCase();
          expect(response.statusMessage).toBe(
            'OK, see redirection instructions in the body of the message');
          let url = `http://localhost/sample?accession=${id}&format=${formatUpperCase}&region=chr1%3A5-401`;
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
          'Range end should be bigger that start');
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

  beforeAll(function() {
    server.listen(socket, () => {
      console.log(`Server listening on socket ${socket}`);
    });
    fse.ensureDirSync(tmpDir);
    options = config.provide(dummy);
  });
  afterAll(function() {
    server.close();
    try { fs.unlinkSync(socket); } catch (e) {}
    fse.removeSync(tmpDir);
  });

  it('data format driven content type', function(done) {
    server.on('request', (request, response) => {
      let c = new RangerController(request, response, {one: "two"}, null, true);
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
