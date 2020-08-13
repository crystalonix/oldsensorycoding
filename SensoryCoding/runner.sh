#! /bin/csh -f

rm sources_list.txt
find ./src/sensoryCoding -name "*.java" > sources_list.txt
javac -cp "./lib/*" @sources_list.txt
cd ./src
java -cp .:"../lib/*" org.junit.runner.JUnitCore sensoryCoding.test.DriverTest >$1 &
