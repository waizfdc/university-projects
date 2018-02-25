#include <boost/graph/use_mpi.hpp>
#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/config.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/strong_components.hpp>

#include <string.h>
#include <vector>


using namespace boost;
using boost::graph::distributed::mpi_process_group;
typedef adjacency_list< vecS, distributedS < mpi_process_group, vecS >, bidirectionalS > Graph;

int main(int argc, char ** argv)
{
	boost::mpi::environment env(argc, argv);
	std::fstream file(argv[1], std::ios::in | std::ios::binary);
	int vertices_count = 0;
	long long edges_count = 0;
	
	// считывание графа
	file.read((char*)(&vertices_count), sizeof(int));
	file.read((char*)(&edges_count), sizeof(long long));
	
	Graph g(vertices_count);
	
	if (process_id(process_group(g)) == 0)
	{
		std::cout << "Graph has " << vertices_count << " vertices" << std::endl;
		std::cout << "Graph has " << edges_count << " edges" << std::endl;
		for(long long i = 0; i < edges_count; i++)
		{
			int src_id = 0, dst_id = 0;
			
			// read i-th edge data
			file.read((char*)(&src_id), sizeof(int));
			file.read((char*)(&dst_id), sizeof(int));

			add_edge(vertex(src_id, g), vertex(dst_id, g), g);
		}
	}
	file.close();


	std::vector<int> local_components_vec(num_vertices(g));
	typedef iterator_property_map<std::vector<int>::iterator, property_map<Graph, vertex_index_t>::type> ComponentMap;
	ComponentMap component(local_components_vec.begin(), get(vertex_index, g));

	// замер времени
	double start = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);

	// запуск подсчёта
	int numComp = strong_components(g, component);


	MPI_Barrier(MPI_COMM_WORLD);
	// подсчёт производительности
	double end = MPI_Wtime();
	if (process_id(process_group(g)) == 0)
	{
		std::cout << num_processes(process_group(g)) << " CPUs; " << "MPI time = " << end - start << " seconds. " 
			<< numComp << " components number; " << "TEPS = " << (1.0 * edges_count)/(end - start) << std::endl;
	}


bool checkWithSeq = false;

// проверка последовательным алгоритмом
if (checkWithSeq)
{
	if ( process_id(process_group(g)) == 0 )
	{
		for (int i = 0; i < vertices_count; ++i)
			get(component, vertex(i, g));
		synchronize(component);

		// Check against the sequential version
		typedef adjacency_list<vecS, vecS, directedS> Graph2;
	  
		file.open(argv[1], std::ios::in | std::ios::binary);
		file.read((char*)(&vertices_count), sizeof(int));
		file.read((char*)(&edges_count), sizeof(long long));

		Graph2 g2(vertices_count);
		
		for(long long i = 0; i < edges_count; i++)
		{
			int src_id = 0, dst_id = 0;
			
			// read i-th edge data
			file.read((char*)(&src_id), sizeof(int));
			file.read((char*)(&dst_id), sizeof(int));
			add_edge(vertex(src_id, g2), vertex(dst_id, g2), g2);
		}		
		file.close();	
		

		std::vector<int> component2(vertices_count);
		int seqNumComp = strong_components(g2, make_iterator_property_map(component2.begin(), get(vertex_index, g2)));

		// проверяем, что результат совпадает
		std::map<int, int> c2c;
		int i;
		for ( i = 0; i < vertices_count; i++ )
			if ( c2c.find( get(component, vertex(i, g)) ) == c2c.end() )
				c2c[get(component, vertex(i, g))] = component2[i];
			else
				if ( c2c[get(component, vertex(i, g))] != component2[i] )
					break;
	  
		if (( i < vertices_count ) || (numComp != seqNumComp))
			std::cout << "Unable to checkWithSeq SCC result...\n";
		else
			std::cout << "Checked! " << seqNumComp << " strong components\n"; 
	}
	else 
	{
		synchronize(component);
	}
}
return 0;
}