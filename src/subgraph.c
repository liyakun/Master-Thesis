#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include "../subgraphIsomorphism.h"
#include "../frequentsubgraph.h"
#include "../graph.h"
#include "../loading.h"

int printHelp(){
#include "./subgraphHelp.help"
	unsigned char* help = executables_subgraphHelp_txt;
	int len = executables_subgraphHelp_txt_len;
	if (help != NULL){
		int i=0;
		for (i=0; i<len; ++i){
			fputc(help[i], stdout);
		}
		return EXIT_SUCCESS;
	}else{
		fprintf(stderr, "Could not read help file\n");
		return EXIT_FAILURE;
	}
}

void subGraphIsomophism(Vector* DB, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	if(DB->size != 2){
		printf("subGraphIsomophism assumes tree in DB[0], and graph in DB[1], and total size of DB is 2!\n");
		printf("You DB has a size %i\n", DB->size);
	}else{
	    clock_t tic = clock();
		subGraphIsoCheck((struct Graph*)(DB->data[0]), (struct Graph*)(DB->data[1]), NULL, gp, sgp);
	    clock_t toc = clock();
	    printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		for(int i=0; i<DB->size; i++){
			dumpGraph(gp, DB->data[i]);
		}
		vector_free(DB);
	}
}

void frequentSubTreeMining(Vector* DB, int support, int order, char miningMethod, char search_method, char subTreeIsomorphic, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	/* pre-processing graphs */
	Vector objects;
	vector_init(&objects);
	graphPreProcessing(DB, &objects, gp, sgp);

	clock_t tic = clock();
	if(miningMethod == 'l'){
		if(subTreeIsomorphic == 'i'){
			apFrequentSubGraphMining(DB, &objects, gp, sgp, support, order, 'f', search_method);
		}else{
			apFrequentSubGraphMining(DB, &objects, gp, sgp, support, order, ' ', search_method);
		}
	}else{
		if(subTreeIsomorphic == 'i'){
			fpFrequentSubGraphMining(DB, &objects, gp, sgp, support, order, 'f', search_method);
		}else{
			fpFrequentSubGraphMining(DB, &objects, gp, sgp, support, order, ' ', search_method);
		}
	}
	clock_t toc = clock();
	printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

	/* post-processing graphs */
	graphPostProcessing(DB, &objects, gp, sgp);
}

int main(int argc, char** argv){

	/* object pool */
	struct ListPool* lp;
	struct VertexPool* vp;
	struct ShallowGraphPool* sgp;
	struct GraphPool* gp;

	/* pointer to the current graph which is returned by the input iterator */
	struct Graph* g = NULL;

	/* user input handling variables */
	int support = 1;
	char miningMethod = 'l';
	char subTreeIsomorphic = 'i'; 
	char search_method = 'p';
	char usedFor = 'm';
	int order = 7;

	/* parse command line arguments */
	int arg;
	const char* validArgs = "hu:m:i:s:v:c:";
	for (arg=getopt(argc, argv, validArgs); arg!= -1; arg=getopt(argc, argv, validArgs)){
		switch(arg){
		case 'h':
			printHelp();
			return EXIT_SUCCESS;
			break;
		case 'u':
			if(strcmp(optarg, "isomorphic") == 0){
				usedFor = 'i';
				break;
			}
			if(strcmp(optarg, "mining") == 0){
				usedFor = 'm';
				break;
			}
			fprintf(stderr, "Unknow used for: %s\n", optarg);
			return EXIT_FAILURE;
			break;
		case 'm':
			if (strcmp(optarg, "levelwise") == 0){
				miningMethod = 'l';
				break;
			}
			if (strcmp(optarg, "frequentpattern") == 0){
				miningMethod = 'f';
				break;
			}
			fprintf(stderr, "Unknow mining method: %s\n", optarg);
			return EXIT_FAILURE;
			break;
		case 'i':
			if(strcmp(optarg, "noniterative") == 0){
				subTreeIsomorphic = 'n';
				break;
			}
			if(strcmp(optarg, "iterative") == 0){
				subTreeIsomorphic = 'i';
				break;
			}
			fprintf(stderr, "Unknown subtree isomorphic method: %s\n", optarg);
			return EXIT_FAILURE;
			break;
		case 's':
			if (sscanf(optarg, "%i", &support) != 1){
				fprintf(stderr, "support value must be number, is: %s\n", optarg);
				return EXIT_FAILURE;
			}
			break;
		case 'v':	
			if (sscanf(optarg, "%i", &order) != 1){
				fprintf(stderr, "order value must be number, is: %s\n", optarg);
				return EXIT_FAILURE;
			}
			break;
		case 'c':
			if (strcmp(optarg, "binary") == 0){
				search_method = 'y';
				break;
			}
			if(strcmp(optarg, "prefix") == 0){
				search_method = 'p';
				break;
			}
			fprintf(stderr, "Unknown candidate elimination search method: %s\n", optarg);
			return EXIT_FAILURE;
			break;
		case '?':
			return EXIT_FAILURE;
			break;
	    }
    }

	/* init object pools */
	lp = createListPool(100000);
	vp = createVertexPool(100000);
	sgp = createShallowGraphPool(100000, lp);
	gp = createGraphPool(1000000, vp, lp);

	/* initialize the stream to read graphs from check if there is a filename present
	 * in the command line arguments, if so, open the file, if not, read from stdin*/
	if (optind < argc){
		char* filename = argv[optind];
		/* if the present filename is not '-' then init a file iterator for that file name*/
		if (strcmp(filename, "-") != 0){
			createFileIterator(filename, gp);
		} else {
			createStdinIterator(gp);
		}
	} else {
		createStdinIterator(gp);
	}

	/* iterate over all graphs in the database */
	printf("\n******************** Load Graph **********************\n");
	Vector DB;
	vector_init(&DB);
	int number_of_graphs = 0;
	while((g = iterateFile())){ //  && number_of_graphs < 529
		if (g->n > 0){
			++number_of_graphs;
			vector_append(&DB, g);
		}
	}

	printf("Number of graphs: %i\n", DB.size);
	if(usedFor == 'i'){
		subGraphIsomophism(&DB, gp, sgp);
	}else{
		frequentSubTreeMining(&DB, support, order, miningMethod, search_method, subTreeIsomorphic, gp, sgp);
	}

	/* global garbage collection */
	destroyFileIterator();
	freeGraphPool(gp);
	freeShallowGraphPool(sgp);
	freeListPool(lp);
	freeVertexPool(vp);

	return EXIT_SUCCESS;
}
