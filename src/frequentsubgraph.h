#ifndef FREQUENTSUBGRAPH_H
#define FREQUENTSUBGRAPH_H
#define VECTOR_INIT_CAPACITY 50
#define DVC_INIT_CAPACITY 20
#define VECTOR_CHARACTERISTICS_INIT_CAPACITY 20
#include "./graph.h"

/* Add vector simulation for dynamic size pointer array
 * referencing http://eddmann.com/posts/implementing-a-dynamic-vector-array-in-c/
 **/
typedef struct vector{
	void** data;
	int capacity;
	int size;
} Vector;

/* Define characteristics structure
 **/
struct CoreLimb{
	struct ShallowGraph* core;
	Vector* limbs; // contain sequential edge-vertex label for core
};

/* Add vector simulation for dynamic size int array
 **/
typedef struct vectorint{
	int* data;
	int capacity;
	int size;
} VectorInt;

/* Add vector of VectorInt
 **/
typedef struct dvectorint{
	int capacity;
	int size;
	VectorInt** data;
} DVectorInt;

/* Add vector simulation for dynamic size characteristics array
 **/
typedef struct vcharacteristics{
	int size;
	int capacity;
	struct IsoCharacteristics* data;
} VectorCharacteristics;

/* Add double vector simulation for dynamic size characteristics array
 **/
typedef struct dcharacteristics{
	int size;
	int capacity;
	VectorCharacteristics** data;
} DVectorCharacteristics;

void dvint_init(DVectorInt*);
void dvint_resize(DVectorInt*, int);
void dvint_realize(DVectorInt*);
void dvint_clean(DVectorInt*);
void dvint_free(DVectorInt*);

void dvc_init(DVectorCharacteristics*);
void dvc_resize(DVectorCharacteristics*, int);
void dvc_realize(DVectorCharacteristics*, int);
void dvc_free(DVectorCharacteristics*, int);

void vint_init(VectorInt*);
void vint_resize(VectorInt*, int);
void vint_append(VectorInt*, int);
int vint_equal(VectorInt*, VectorInt*);
VectorInt* vint_copy(VectorInt*);
int vint_distinct(VectorInt*);
void vint_free(VectorInt*);

void vc_init(VectorCharacteristics*);
void vc_resize(VectorCharacteristics*, int);
void vc_append(VectorCharacteristics*, struct IsoCharacteristics*);
int vc_contain(VectorCharacteristics*, struct IsoCharacteristics*);
void vc_free(VectorCharacteristics*);

void vector_init(Vector*);
void vector_resize(Vector*, int);
void vector_append(Vector*, void*);
void vector_free(Vector*);
void vector_copy(Vector*, Vector*);

/* Create a graph with single vertex, and labeled with "label"
 **/
struct Graph* createSingleVertexGraph(struct GraphPool* gp, char* label);

/* Pre process graphs used for frequent tree mining
 * 0. get all the 2 connected components in each graph
 * 1. get all joined 2 connected components
 * 2. get all locally generated spanning trees
 * 3. get all post orders in each locally generated spannning trees
 **/
void graphPreProcessing(Vector* graphDB, Vector* graphObjects, struct GraphPool* gp, struct ShallowGraphPool* sgp);

/* Post process the list of generated objects used in frequent subtree mining
 * this part frees all objects corresponding to the generated objects in graphPreProcessing()
 **/
void graphPostProcessing(Vector* graphDB, Vector* graphObjects, struct GraphPool* gp, struct ShallowGraphPool* sgp);

/* Implement Algorithm 1 from paper ilp2014long-cameraReady_v2.pdf, with apriori fashion
 * 1. flag == 'f', means use fastSubGraphIsoCheck, algorithm from main.pdf
 * 2. flag == ' ', means use algorithm subGraphIsoCheck, algorithm 1 Subgraph Isomorphism from a Tree into a Connected Graph, ipl2015-submitted-draft.pdf
 **/
void apFrequentSubGraphMining(Vector* DB, Vector* objects, struct GraphPool* gp, struct ShallowGraphPool* sgp, int t, int order, char flag, char search);

/* Implement Algorithm 1 from paper ilp2014long-cameraReady_v2.pdf, with fp growth,
 * 1. flag == 'f', means use fastSubGraphIsoCheck, algorithm from main.pdf
 * 2. flag == ' ', means use algorithm subGraphIsoCheck, algorithm 1 Subgraph Isomorphism from a Tree into a Connected Graph, ipl2015-submitted-draft.pdf
 **/
void fpFrequentSubGraphMining(Vector* GraphDB, Vector* objects, struct GraphPool* gp, struct ShallowGraphPool* sgp, int t, int order, char flag, char search);
#endif
