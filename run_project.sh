# rm BogdanMain
# rm Main1
# rm Main2
rm MainBFV
# rm MainCKKS

# g++ -g -O2 -std=c++11 -pthread -march=native BogdanMain.cpp -o BogdanMain -lntl -lgmp -lm
# g++ -g -O2 -std=c++11 -pthread -march=native Main1.cpp -o Main1 -lntl -lgmp -lm
# g++ -g -O2 -std=c++11 -pthread -march=native Main2.cpp -o Main2 -lntl -lgmp -lm
g++ -g -O2 -std=c++11 -pthread -march=native MainBFV.cpp -o MainBFV -lntl -lgmp -lm
# g++ -g -O2 -std=c++11 -pthread -march=native MainCKKS.cpp -o MainCKKS -lntl -lgmp -lm

# ./BogdanMain
# ./Main1
# ./Main2
./MainBFV
# ./MainCKKS


