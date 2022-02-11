> $1 # clear file

for i in {1..10}
do
    echo "${i} / 10"
    bazel-bin/sph/sph /home/jay/Dropbox/University/Y4/PHYM004/sph/config.txt >> $1
done
