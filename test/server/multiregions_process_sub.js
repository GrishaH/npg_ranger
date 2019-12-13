/*globals describe, it, expect*/

"use strict";


const multiRegions = require('../../lib/server/Ztest_multi_regions.js');

describe('Test valid multiregion string', function() {
  
  it("Basic single view and merge", function() {
    let testingString = multiRegions.createCommand(["testdata1.bam"], []);
    console.log(testingString);
    expect(testingString).toBe("samtools merge -O SAM - <(samtools view -h testdata1.bam )");
  });

  it("Basic single view and merge with given region", function() {
    let testingString = multiRegions.createCommand(["testdata1.bam"], ["phix:300-500"]);
    console.log(testingString);
    expect(testingString).toBe("samtools merge -O SAM - <(samtools view -h testdata1.bam phix:300-500 )");
  });

  it("Basic single view and merge with several regions", function() {
    let testingString = multiRegions.createCommand(["testdata1.bam"], ["phix:300-400", "phix:600-1200"]);
    console.log(testingString);
    expect(testingString).toBe("samtools merge -O SAM - <(samtools view -h testdata1.bam phix:300-400 phix:600-1200 )");
  });

  it("View three files and merge", function() {
    let testingString = multiRegions.createCommand(["testdata1.bam", "testdata2.bam", "testdata3.bam"], []);
    console.log(testingString);
    expect(testingString).toBe("samtools merge -O SAM - <(samtools view -h testdata1.bam ) <(samtools view -h testdata2.bam ) <(samtools view -h testdata3.bam )");
  });

  it("View three files and merge with given region", function() {
    let testingString = multiRegions.createCommand(["testdata1.bam", "testdata2.bam", "testdata3.bam"], ["phix:300-500"]);
    console.log(testingString);
    expect(testingString).toBe("samtools merge -O SAM - <(samtools view -h testdata1.bam phix:300-500 ) <(samtools view -h testdata2.bam phix:300-500 ) <(samtools view -h testdata3.bam phix:300-500 )");
  }); 

  it("View three files and merge with two given regions", function() {
    let testingString = multiRegions.createCommand(["testdata1.bam", "testdata2.bam", "testdata3.bam"], ["phix:300-500", "phix:600-1200"]);
    console.log(testingString);
    expect(testingString).toBe("samtools merge -O SAM - <(samtools view -h testdata1.bam phix:300-500 phix:600-1200 ) <(samtools view -h testdata2.bam phix:300-500 phix:600-1200 ) <(samtools view -h testdata3.bam phix:300-500 phix:600-1200 )");
  }); 

  /*
  it("No inputs given", function() { // ?
    let testingString = multiRegions.createCommand([""]);
    console.log(testingString);
    expect(testingString).toBe("samtools merge -O SAM - <(samtools view -h testdata1.bam )");
  }); 
  */

});