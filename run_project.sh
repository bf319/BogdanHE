# rm BogdanMain
# rm Main1
# rm Main2
# g++ -g -O2 -std=c++11 -pthread -march=native BogdanMain.cpp -o BogdanMain -lntl -lgmp -lm
# g++ -g -O2 -std=c++11 -pthread -march=native Main1.cpp -o Main1 -lntl -lgmp -lm
# g++ -g -O2 -std=c++11 -pthread -march=native Main2.cpp -o Main2 -lntl -lgmp -lm
# ./BogdanMain
# ./Main1
# ./Main2

rm MainFinal
g++ -g -O2 -std=c++11 -pthread -march=native MainFinal.cpp -o MainFinal -lntl -lgmp -lm
./MainFinal