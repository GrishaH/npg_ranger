"use strict";

const fs = require("fs");
const spawn = require('child_process').spawn;

let testFilenameArray = ["testdata1.bam", "testdata2.bam"];
let regionArray = ["phix:200-300", "phix:400-1000"];


let createCommand = (fileNames, regionArray) => {
  let finalString =  "samtools merge -O SAM -";
  let regionString = " ";
  // Form the string used for regions
  for (let region of regionArray) {
    regionString += region + " ";
  }
  for (let filename of fileNames) {
    finalString += " <(samtools view -h " + filename + regionString  + ")";
  }
  // console.log(finalString);
  return finalString;
};

var multiregionMerge = (fileNames, regionArray = []) => {
  // Check region array
  if (!Array.isArray(regionArray)) {
    console.log("Not an array!");
    regionArray = [];
    // Either throw error or fix array?
  }
  let inputShellCommand = createCommand(fileNames, regionArray);
  // let inputShellCommand = "diff <(cat data1.txt) <(cat data2.txt)";

  let writeStream = fs.createWriteStream("outputFile.txt");

  let child = spawn('/bin/bash', ['-c', inputShellCommand], { stdio: 'pipe' }); // pipe or inherit
  child.stdout.pipe(writeStream);
  // child.stderr.pipe(writeStream);

  child.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    throw data;
  });

  child.on('close', (code) => {
    console.log(`child process exited with code ${code}`);
  });

  // child.stdout.on('data', (data) => {
  //   console.log(`stdout: ${data}`);
  // });
  /*
  child.on('exit', (code, signal) => {
    if (code) {
      console.error('Child exited with code', code)
    } else if (signal) {
      console.error('Child was killed with signal', signal);
    } else {
      console.log('Child exited okay');
    }
  });

  child.stdout.on('data', (data) => {
    console.log('-Writing data-');
    writeStream.write(data); // pipe instead of data?
    // console.log(`stdout:\n${data}`);
  });
  */
};

multiregionMerge(testFilenameArray, regionArray);

module.exports = {
  multiregionMerge: multiregionMerge,
  createCommand: createCommand
};


/*
// use TESTJSON2.json and TESTJSON3.json
const spawn = require('child_process').spawn;
const { exec } = require('child_process');

let spawnOptions = {
  stdio: []
};

let test1 = spawn('cat', ['TESTJSON2.json']);

test1.stdout.on('data', (data) => {
  console.log(`stdout: ${data}`);
});


let p = execSync('"diff <(cat TESTJSON2.json) <(cat TESTJSON3.json)"',
{shell: '/bin/bash'}, (error, stdout, stderr) => {
  if (error) {
    console.error(`error: ${error}`);
    return;
  }
  console.log(`stdout: ${stdout}`);
  console.error(`stderr: ${stderr}`);
});
*/

/*
let readFile = (name) => {
  return new Promise((resolve, reject) => {
    fs.readFile(name, function(err, data){
      if (err) {
        reject(err);
      } else {
        resolve(data);
      }
    });
  });
};
*/



/*
let p = exec("diff <(cat a.txt) <(cat b.txt)", {shell: '/bin/bash'},
(error, stdout, stderr) => {
  if (error) {
    console.error(`exec error: ${error}`);
    return;
  }
  console.log(`stdout: ${stdout}`);
  console.error(`stderr: ${stderr}`);
});
*/

/*
Promise.all([readFile("read_data_1"), readFile("read_data_2")])
.catch(err => {
  console.log(err);
})
.then(data => {
  let p = exec("diff " + data[0] + data[1],
               (error, stdout, stderr) => {
    if (error) {
      console.error(`exec error: ${error}`);
      return;
    }
    console.log(`stdout: ${stdout}`);
    console.error(`stderr: ${stderr}`);
  });
  console.log(data[0]);
  console.log(data[1]);
});
*/




/*
console.log("Before execution");
fs.readFile("read_data_1", (err, data) => {
  console.log("Reading file");
  if (err) {
    console.log("Error: " + err);
  }
  console.log("Avoided err");
  console.log(err);
  console.log(data.toString());
});
*/