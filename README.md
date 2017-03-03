#WES QC Hail NOTES


###Andrea Ganna's steps for running QC and analysis in WGS data:

[Analysis Plan](https://storage.googleapis.com/wgspd-wgs-v1-noshared/Analysis_plan_METSIM.md.html)

[WGSPD] (https://storage.googleapis.com/wgspd-wgs-v1/wgspd_guide_shared.html)



###Hail VEP from Cotton:

```
#!/bin/bash

export PATH=/usr/java/jdk1.8.0_60/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin
export JAVA_HOME=/usr/java/jdk1.8.0_60

. /psych/genetics_data/working/cseed/hail-inst/etc/jar.sh

export SPARK_CLASSPATH=$JAR

# 8g/core overhead
/usr/bin/spark-submit --conf spark.yarn.max.executor.failures=500 --conf spark.task.maxFailures=50 --conf spark.yarn.executor.memoryOverhead=40960 --driver-memory 20g --executor-memory 10g --executor-cores 5 --class org.broadinstitute.hail.driver.Main $JAR --master yarn-client "$@"
```
