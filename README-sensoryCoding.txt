# Project Title
SensoryCoding: The objective of this project is to build an end-to-end coding/decoding framework that transforms continuous time signals into spike trains and finally reconstructs the signal back. As of now this project only deals with one dimensional
Signal.

## Getting Started

First you need to install the project in your local machine. For that you need to install git in your local machine. Once git is installed run the following command to get a copy of the project in your local machine:
git clone https://bitbucket.org/crystalonix/oldsensorycoding.git

### Prerequisites

What things you need to install the software and how to install them
For running the code you need installed jdk version 8 or above.

### Installing

2. First download the following jar files:
ejml-cdense-0.32-sources.jar; 		ejml-zdense-0.32-sources.jar;
ejml-cdense-0.32.jar;			ejml-zdense-0.32.jar;
ejml-core-0.32-sources.jar;		hamcrest-core-1.3.jar;
ejml-core-0.32.jar;			jcommon-1.0.23.jar;
ejml-ddense-0.32-sources.jar;		jfreechart-1.0.19-experimental.jar;
ejml-ddense-0.32.jar;			jfreechart-1.0.19-swt.jar;
ejml-dsparse-0.32-sources.jar;		jfreechart-1.0.19.jar;
ejml-dsparse-0.32.jar;			jfreesvg-2.0.jar;
ejml-experimental-0.32-sources.jar;	junit-4.11.jar;
ejml-experimental-0.32.jar;		orsoncharts-1.4-eval-nofx.jar;
ejml-fdense-0.32-sources.jar;		orsonpdf-1.6-eval.jar;
ejml-fdense-0.32.jar;			servlet.jar;
ejml-simple-0.32-sources.jar;		swtgraphics2d.jar;
ejml-simple-0.32.jar;

(Most of the jars correspond to ejml and jfreechart library dependency, and some few extra jars are needed)

2. Once all the jars are downloaded put them in the lib folder of the project;

3. Go to the root directory of the project and run the following commands to compile all the .java files:

find ./src/sensoryCoding -name "*.java" > sources_list.txt
javac -cp "./lib/*" @sources_list.txt

## Running the tests

All tests are included inside the test folder. Before running any test please make the necessary configuration changes. Configurations can be changed explicitly by rewriting the entries in ConfigurationParameters.java file or some configurations are embedded in the test class method itself (e.g. location of the extracted audio file). Samples of extracted audio files are present inside the soundData folder. For running all the tests of a specific test class do the following:
1. cd ./src
2. java -Xmx64G -cp .:"../lib/*" org.junit.runner.JUnitCore sensoryCoding.test.TESTCLASSNAME > LOGFILENAME

Please Note: Some of the operations are memory intensive and hence it recommended to allocate maximum possible java heap space (here we kept -Xmx64G).  

### Break down into end to end tests



We used git for version control.
