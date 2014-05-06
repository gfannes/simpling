.PHONY: test

CPP_OPTIONS := -std=c++11

test: cpp/test/MetropolisHastings_tests.cpp
	g++ $(CPP_OPTIONS) cpp/test/MetropolisHastings_tests.cpp -o MetropolisHastings_tests.exe -Icpp
	./MetropolisHastings_tests.exe
