### Some common tasks in bioinformatics

#### Homolog search

- blat

```{bash}
#!/bin/bash
if [ $# -ne 2 ];then
echo "Usage: make-ooc.sh db ooc"
exit 
fi
db=$1
ooc=$2
blat ${db} /dev/null -q=dnax -t=dnax -stepSize=11  -makeOoc=${ooc} /dev/null
```



#### Remove rRNA from sequencing reads
