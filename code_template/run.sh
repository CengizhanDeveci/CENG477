make all

for FILE in ../inputs/*;
do
    echo Starts $FILE;
    start_time=$(date +%s);
    ./raytracer $FILE;
    end_time=$(date +%s);
    elapsed=$(( end_time - start_time ));
    echo $FILE;
    echo $elapsed;
done