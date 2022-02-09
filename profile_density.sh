> $1 # clear file

for i in {1..3}
do
    bazel-bin/sph/sph /home/jay/Dropbox/University/Y4/PHYM004/sph/config.txt >> $1
done
