all: MergeMap.exe

MergeMap.exe: main.o single_population_1.o single_population_2.o map_integration.o condense_markers_into_bins.o fragment.o ConsensusLG.o SCC_Drawer.o shortest_path.o heuristic_solver.o
	g++ -o consensus_map.exe  main.o single_population_1.o single_population_2.o map_integration.o condense_markers_into_bins.o fragment.o ConsensusLG.o SCC_Drawer.o shortest_path.o heuristic_solver.o

main.o: main.cpp constants.h
	g++ -o main.o -c -g main.cpp -I/home/yonghui/boost/boost_1_34_1/include/boost-1_34_1

single_population_1.o : single_population.cpp single_population.h 
	g++ -o single_population_1.o -c -g single_population.cpp  -I/home/yonghui/boost/boost_1_34_1/include/boost-1_34_1

single_population_2.o : single_population_linearize_dag.cpp single_population.h 
	g++ -o single_population_2.o -c -g single_population_linearize_dag.cpp -I/home/yonghui/boost/boost_1_34_1/include/boost-1_34_1


map_integration.o : map_integration.cpp map_integration.h
	g++ -o map_integration.o -c -g map_integration.cpp -I/home/yonghui/boost/boost_1_34_1/include/boost-1_34_1 
	
condense_markers_into_bins.o : condense_markers_into_bins.cpp condense_markers_into_bins.h
	g++ -o condense_markers_into_bins.o -c -g condense_markers_into_bins.cpp -I/home/yonghui/boost/boost_1_34_1/include/boost-1_34_1

fragment.o : fragment.cpp fragment.h
	g++ -o fragment.o -c -g fragment.cpp -I/home/yonghui/boost/boost_1_34_1/include/boost-1_34_1
	
ConsensusLG.o : ConsensusLG.cpp ConsensusLG.h
	g++ -o ConsensusLG.o -c -g ConsensusLG.cpp -I/home/yonghui/boost/boost_1_34_1/include/boost-1_34_1	

SCC_Drawer.o : SCC_Drawer.cpp SCC_Drawer.h
	g++ -o SCC_Drawer.o -c -g SCC_Drawer.cpp -I/home/yonghui/boost/boost_1_34_1/include/boost-1_34_1

shortest_path.o : shortest_path.cpp shortest_path.h
	g++ -o shortest_path.o -c -g shortest_path.cpp -I/home/yonghui/boost/boost_1_34_1/include/boost-1_34_1

heuristic_solver.o : heuristic_solver.cpp heuristic_solver.h
	g++ -o heuristic_solver.o -c -g heuristic_solver.cpp 
	
clean:
	rm -rf *.o; rm -rf *.exe; rm -rf *~; rm -rf *.out
