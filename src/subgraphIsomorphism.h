#ifndef SUBGRAPHISOMORPHISM_H
#define SUBGRAPHISOMORPHISM_H
#include "./graph.h"
#include "frequentsubgraph.h"

/* Define characteristics structure
 **/
struct IsoCharacteristics{
	int start;	// y vertex number
	int end;	// u vertex number
	struct ShallowGraph* t; // pointer of local generated spanning tree
	int w;		// w vertex in graph G
};

struct GraphObjects{
	struct Graph *g;
	Vector *joi_bi_components;
	Vector *spanning_trees;
	Vector *post_orders;
	Vector *sub_graphs;
	struct ShallowGraph *b_components;
	struct GraphObjects *next;
};

struct CandidateObjects{
	struct Graph* h;
	VectorInt *graphs;
	Vector *characteristics;
	//char* canonicalString;
	struct ShallowGraph* canonicalString;
	struct CandidateObjects* parent;
	int number;
};

/* get the vertex number by id
 **/
int getNumById(struct Graph* g, int v_id);

int subtreeIsomorphic(VectorCharacteristics* CC, struct Graph* H, struct Graph* G, struct GraphObjects* object);

/* Implement algorithm 1 in paper "Subtree Isomorphism in Graphs with Locally Polynomial Spanning Trees"
	by Pascal Welke, Tamás Horváth, Stefan Wrobel
 **/
int subGraphIsoCheck(struct Graph* H, struct Graph* G, struct GraphObjects*, struct GraphPool* gp, struct ShallowGraphPool* sgp);

/* Do subgraph isomorphic check from H to G, by using information from (H-e)
 **/
int fastSubgraphIsoCheck(struct Graph* H, struct Graph* G, struct GraphPool* gp, struct ShallowGraphPool* sgp);

/* Implement Algorithm 2 from main.pdf
 **/
int fastSubGraphIsoCheck(struct Graph* H, struct Graph* G, struct GraphObjects*, VectorCharacteristics* CC_pre, VectorCharacteristics* CC_current, int a, int b);

/* Join 2-connected components with common cut vertex */
void joinBiComponents(Vector* joi_bi_components, struct ShallowGraph* bComponents, struct Graph* G);

/* Get the post order of spanning tree t, rooted vertex v
 **/
VectorInt* getSpanningTreePostOrder(struct ShallowGraph* t, struct GraphPool* gp, struct Vertex* v);

/* Initialize characteristics cha
 **/
void initCharacteristic(struct IsoCharacteristics* cha);

/* Create a new graph g by keeping vertices of graph G in vertices list vl
 **/
struct Graph* createSimpleGraph(struct GraphPool* gp, struct Graph* H, VectorInt* vl);

/* Create a new graph g by removing vertex 'vertex_remove' from graph G
 **/
struct Graph* createSGraph(struct GraphPool* gp, struct VertexPool* vp, struct Graph* G, int vertex_remove);

/* Pre processing a single graph
 * 0. get all the 2 connected components in each graph
 * 1. get all joined 2 connected components
 * 2. get all locally generated spanning trees
 * 3. get all post orders in each locally generated spanning trees
 */
struct GraphObjects* singleGraphPreProcessing(struct Graph* g, struct ShallowGraph* bComponents, Vector* joi_bi_components, Vector* spanning_tree_list, Vector* sub_graphs, Vector* post_orders, struct GraphPool* gp, struct ShallowGraphPool* sgp);

/* Post process the list of generated objects used in graph G, used in subtree isomorphic
 * this part frees all objects corresponding to the generated objects in singleGraphPreProcessing()
 */
void singleGraphPostProcessing(struct GraphObjects* object, struct GraphPool* gp, struct ShallowGraphPool* sgp);

/* Print the characteristics in CC
 **/
void printCC(VectorCharacteristics* CC);
#endif
