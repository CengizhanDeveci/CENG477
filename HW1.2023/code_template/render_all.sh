make clean
make
for file in ../inputs/*
do
    echo $file
    ./raytracer $file
done
