#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include "./subgraphIsomorphism.h"
#include "./frequentsubgraph.h"
#include "./graphPrinting.h"
#include "./listComponents.h"
#include "./listSpanningTrees.h"
#include "./cs_Tree.h"
#include "./cs_Parsing.h"
#include "./cs_Compare.h"
#include "./searchTree.h"

/* Initialize vector v with initial size
 **/
void vector_init(Vector* v)
{
	v->capacity = VECTOR_INIT_CAPACITY;
	v->size = 0;
	v->data = malloc(sizeof(void*) * v->capacity);
}

/* Resize vector v with specified capacity
 **/
void vector_resize(Vector* v, int capacity)
{
    // printf("vector_resize: %d to %d\n", v->capacity, capacity);
	void** items = realloc(v->data, sizeof(void*) * capacity);
	if (items) {
		v->data = items;
		v->capacity = capacity;
	}
}

/* Append item to current vector v, if size exceeds capacity, double capacity
 **/
void vector_append(Vector* v, void* item)
{
	if (v->capacity == v->size){
		vector_resize(v, v->capacity * 2);
	}
	v->data[v->size++] = item;
}

/* Copy elements from v2 to v1
 **/
void vector_copy(Vector* v1, Vector* v2){
	for(int i=0; i<v2->size; i++){
		vector_append(v1, v2->data[i]);
	}
}

void vector_insert(Vector* v, int index, void* p){
	if(index >= 0 && index <= v->size){
		if(v->capacity == v->size){
			vector_resize(v, v->capacity * 2);
		}

		void* tmp = v->data[index];
		v->data[index] = p;

		for(int i=index+1; i<=v->size; i++){
			void* temp = v->data[i];
			v->data[i] = tmp;
			tmp = temp;
		}

		++v->size;
	}else{
		printf("Error to insert to index: %i\n", index);
	}
}

/* Free vector v
 **/
void vector_free(Vector* v)
{
	free(v->data);
}

/* Initialize vector of integer type
 **/
void vint_init(VectorInt* vint){
	vint->capacity = VECTOR_INIT_CAPACITY;
	vint->size = 0;
	vint->data = (int*)malloc(sizeof(int) * vint->capacity);
}

/* Resize a vector of integer type
 **/
void vint_resize(VectorInt* vint, int capacity){
	int* items = (int*)realloc(vint->data, sizeof(int) * capacity);
	if(items){
		vint->data = items;
		vint->capacity = capacity;
	}
}

/* Append integer item to current integer vector
 **/
void vint_append(VectorInt* vint, int item){
	if(vint->capacity == vint->size){
		vint_resize(vint, vint->capacity * 2 );
	}
	vint->data[vint->size++] = item;
}

void vint_insert(VectorInt* vint, int index, int item){
	if(index >= 0 && index <= vint->size){
		if(vint->capacity == vint->size){
			vint_resize(vint, vint->capacity * 2);
		}

		int tmp = vint->data[index];
		vint->data[index] = item;

		for(int i=index+1; i<=vint->size; i++){
			int temp = vint->data[i];
			vint->data[i] = tmp;
			tmp = temp;
		}

		++vint->size;
	}else{
		printf("Error to insert to index: %i\n", index);
	}
}

/* Check whether vint1 and vint2 is equal sequentially
 **/
int vint_equal(VectorInt* vint1, VectorInt* vint2){
	if(vint1->size != vint2->size){
		return 0;
	}else{
		for(int i=0; i<vint1->size; i++){
			if(vint1->data[i] != vint2->data[i]){
				return 0;
			}
		}
	}
	return 1;
}

/* Get a copy of vint, free using: vector_free(vintcopy); free(vintcopy)
 **/
VectorInt* vint_copy(VectorInt* vint){
	VectorInt* vintcopy = malloc(sizeof(VectorInt));
	vint_init(vintcopy);
	for(int i=0; i<vint->size; i++){
		vint_append(vintcopy, vint->data[i]);
	}
	return vintcopy;
}

/* Get the distinct value in VectorInt to first NewLengh position
 **/
int vint_distinct(VectorInt* vint){
	int i, j, NewLength = 1;
	for(i=1; i< vint->size; i++){
			for(j=0; j< NewLength ; j++)
			{
				if(vint->data[i] == vint->data[j])
					break;
			}
			if (j==NewLength)
				vint->data[NewLength++] = vint->data[i];
	}
	return NewLength;
}

int vint_contain(VectorInt* vint, int value){
	for(int i=0; i<vint->size; i++){
		if(vint->data[i] == value){
			return i;
		}
	}
	return -1;
}

/* Free vint vector
 **/
void vint_free(VectorInt* vint){
	free(vint->data);
}

/* Initialize vector of VectorInt
 **/
void dvint_init(DVectorInt* dvint){
	dvint->capacity = DVC_INIT_CAPACITY;
	dvint->size = 0;
	dvint->data = (VectorInt**)malloc(sizeof(VectorInt*) * dvint->capacity);
}

/* Resize a vector of VectorInt
 **/
void dvint_resize(DVectorInt* dvint, int capacity){
	VectorInt** items = (VectorInt**)realloc(dvint->data, sizeof(VectorInt*) * capacity);
	if(items){
		dvint->data = items;
		dvint->capacity = capacity;
	}
}

/* Realize VectorInt in DVectorInt
 **/
void dvint_realize(DVectorInt* dvint){
	if(dvint->capacity == dvint->size){
		dvint_resize(dvint, dvint->capacity*2);
	}
	dvint->data[dvint->size] = (VectorInt*)malloc(sizeof(VectorInt));
	vint_init(dvint->data[dvint->size++]);
}

/* Clean last added VectorInt
 **/
void dvint_clean(DVectorInt* dvint){
	vint_free(dvint->data[dvint->size-1]);
	free(dvint->data[dvint->size-1]);
	--(dvint->size);
}

/* Free a vector of VectorInt
 **/
void dvint_free(DVectorInt* dvint){
	for(int i=0; i<dvint->size; i++){
		vint_free(dvint->data[i]);
		free(dvint->data[i]);
	}
	free(dvint->data);
}

/* Initialized a vector of characteristics
 **/
void vc_init(VectorCharacteristics* vc){
	vc->capacity = VECTOR_CHARACTERISTICS_INIT_CAPACITY;
	vc->size = 0;
	vc->data = (struct IsoCharacteristics*)malloc(vc->capacity * sizeof(struct IsoCharacteristics));
}

/* Resize vector of characteristics
 **/
void vc_resize(VectorCharacteristics* vc, int capacity){
	struct IsoCharacteristics* items = (struct IsoCharacteristics*)realloc(vc->data, sizeof(struct IsoCharacteristics) * capacity);
	if (items) {
		vc->data = items;
		vc->capacity = capacity;
	}

}

/* Append characteristics item to current characteristics vector
 **/
void vc_append(VectorCharacteristics* vc, struct IsoCharacteristics* item){
	if(vc->capacity == vc->size){
		vc_resize(vc, vc->capacity * 2);
	}
	vc->data[vc->size].start = item->start;
	vc->data[vc->size].end = item->end;
	vc->data[vc->size].t = item->t;
	vc->data[vc->size].w = item->w;
	++(vc->size);
}

/* Check whether characteristic item is in vector vc
 **/
int vc_contain(VectorCharacteristics* vc, struct IsoCharacteristics* item){
	for(int i=0; i<vc->size; i++){
		if(vc->data[i].start == item->start && vc->data[i].end == item->end && vc->data[i].t == item->t && vc->data[i].w == item->w){
			return 1;
		}
	}
	return 0;
}

/* Free characteristics vector vc
 **/
void vc_free(VectorCharacteristics* vc){
	free(vc->data);
}

/* Initialized a double VectorCharacteristics structure
 **/
void dvc_init(DVectorCharacteristics* dvc){
	dvc->capacity = DVC_INIT_CAPACITY;
	dvc->size = 0;
	dvc->data = malloc(sizeof(VectorCharacteristics*) * dvc->capacity);
}

/* Resize double VectorCharacteristics
 **/
void dvc_resize(DVectorCharacteristics* dvc, int capacity){
	VectorCharacteristics** items = realloc(dvc->data, sizeof(VectorCharacteristics*) * capacity);
	if (items) {
		dvc->data = items;
		dvc->capacity = capacity;
	}
}

/* Realize one element of double VectorCharacteristics
 **/
void dvc_realize(DVectorCharacteristics* dvc, int number_of_graphs){
	if(dvc->capacity == dvc->size){
		dvc_resize(dvc, dvc->capacity * 2);
	}
	dvc->data[dvc->size] = malloc(number_of_graphs * sizeof(VectorCharacteristics));
	for(int j=0; j<number_of_graphs; j++){
		vc_init(&(dvc->data[dvc->size][j]));
	}
	++(dvc->size);
}

/* Clean last realized element of double VectorCharacteristics
 **/
void dvc_clean(DVectorCharacteristics* dvc, int number_of_graphs){
	for(int j=0; j<number_of_graphs; j++){
		vc_free(&(dvc->data[dvc->size-1][j]));
	}
	free(dvc->data[dvc->size-1]);
	--(dvc->size);
}

/* Free a double VectorCharacteristics
 **/
void dvc_free(DVectorCharacteristics* dvc, int number_of_graphs){
	for(int i=0; i<dvc->size; i++){
		for(int j=0; j<number_of_graphs; j++){
			vc_free(&(dvc->data[i][j]));
		}
		free(dvc->data[i]);
	}
	free(dvc->data);
}

/* Pre-process graphs used for frequent tree mining
** by call singleGraphPostProcessing()
**/
void graphPreProcessing(Vector* graphDB, Vector* graphObjects, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	for(int i=0; i<graphDB->size; i++){
		struct ShallowGraph* bComponents_ = NULL;
	    Vector* join = (Vector*)malloc(sizeof(Vector));
	    Vector* spanning_trees = (Vector*)malloc(sizeof(Vector));
    	Vector* post_orders = (Vector*)malloc(sizeof(Vector));
    	Vector* sub_graphs = (Vector*)malloc(sizeof(Vector));
		vector_init(join);
		vector_init(spanning_trees);
		vector_init(sub_graphs);
		vector_init(post_orders);
		struct GraphObjects* obj = singleGraphPreProcessing(graphDB->data[i], bComponents_, join, spanning_trees, sub_graphs, post_orders, gp, sgp);
		vector_append(graphObjects, obj);
	}
}

/* Post process the list of generated objects used in frequent subtree mining
 * this part frees all objects corresponding to the generated objects in graphPreProcessing()
 **/
void graphPostProcessing(Vector* graphDB, Vector* graphObjects, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	// for all the graphs, and graphPreProcessing generated objects
	for(int i=0; i<graphDB->size; i++){
		// free objects in each graph
		singleGraphPostProcessing(graphObjects->data[i], gp, sgp);
		// free graph
		dumpGraph(gp, graphDB->data[i]);
	}
	vector_free(graphObjects);
	/*
	vector_free(spanning_trees_list);
	vector_free(all_post_orders);
	vector_free(all_sub_graphs);
	vector_free(bComponents);
	vector_free(joi_bi_components);*/
	vector_free(graphDB);
}

/* Get a copy of vertex v by keeping label and id of vertex v
 **/
struct Vertex* copyVertex(int n, struct Vertex* v, struct VertexPool* vp){
	struct Vertex *new = getVertex(vp);
	new->number = n;
	new->id = n;
	new->label = v->label;
	return new;
}

/* Add copy of vertex end to the current graph, and add one edge between start and the copied vertex
 * The graph g is changed by this function, with one more vertex, and edges added
 **/
struct Graph* addSingleVertex(struct Graph* g, struct Vertex* start, struct Vertex* end, struct GraphPool* gp, char* label){
	struct Vertex* copy = copyVertex(g->n, end, gp->vertexPool);
	struct Vertex** vertices_original = g->vertices;
	// redirect vertices of graph g to new position
	g->vertices = malloc((g->n + 1) * sizeof(struct Vertex*));
	// copy all the original vertex pointers
	for(int i=0; i<g->n; i++){
		g->vertices[i] = vertices_original[i];
	}
	// add extended vertex pointer
	g->vertices[g->n] = copy;
	// free original vertices pointers of graph g
	free(vertices_original);
	// add edges between start and end vertex, with label, the neighborhood information is also added here, with ++(g->m)
	addEdgeBetweenVertices(start->number, copy->number, label, g, gp);
	// add the number of vertices of graph g
	++(g->n);
	return g;
}

/* Check whether vertex v is in vertex array vertices
 **/
int vertexUnique(struct Vertex* v, Vector* vertices){
	for(int i=0; i<vertices->size; i++){
		if(strcmp(v->label, ((struct Vertex*)vertices->data[i])->label) == 0){	// label is the same
			return i;	// not unique
		}
	}
	return -1;	// label is unique, not in current vertex set
}

/* Check whether edge e is in edge array edges
 **/
int edgeUnique(struct VertexList* e, Vector* edges){
	for(int i=0; i<edges->size; i++){
		struct VertexList* edge = edges->data[i];
		if(strcmp(e->label, edge->label) == 0){
			if(((strcmp(e->startPoint->label, edge->startPoint->label) == 0) && (strcmp(e->endPoint->label, edge->endPoint->label) == 0))
			|| ((strcmp(e->startPoint->label, edge->endPoint->label) == 0) && (strcmp(e->endPoint->label, edge->startPoint->label) == 0))){
				return i; // not unique
			}
		}
	}
	return -1;
}

/* Get vertices from graph database DB with duplicate label vertices removed
 **/
void generateDiffVertices(Vector* GraphDB, Vector* vertices){
	// get different vertices label from the graph DB
	for(int i=0; i<GraphDB->size; i++){
		struct Graph* g = (struct Graph*)GraphDB->data[i];
		for(int j=0; j<g->n; j++){
			struct Vertex* v = g->vertices[j];
			if(vertexUnique(v, vertices)){
				vector_append(vertices, v);
			}
		}
	}
}

void scanVerticesSingleGraph(struct Graph* g, Vector* vertices, DVectorInt* vertices_graphs){
	for(int i=0; i<g->n; i++){
		struct Vertex* v = g->vertices[i];
		int vertex_index = vertexUnique(v, vertices);
		if(vertex_index == -1){ // if current vertices vector did not contain this vertex, then add it
			dvint_realize(vertices_graphs);
			vint_append(vertices_graphs->data[vertices->size], g->number);
			vector_append(vertices, v);
		}else{
			if(vint_contain(vertices_graphs->data[vertex_index], g->number) == -1){
				vint_append(vertices_graphs->data[vertex_index], g->number);
			}
		}
	}
}

void generateVertices(Vector* graphDB, Vector* vertices, DVectorInt* graphsp){
	struct Graph* g = NULL;
	for(int i=0; i<graphDB->size; i++){
		g = graphDB->data[i];
		scanVerticesSingleGraph(g, vertices, graphsp);
	}
}

void scanEdgesSingleGraph(struct Graph* g, Vector* edges, DVectorInt* edges_graphs, int index){
	struct VertexList* e = NULL;
	for(int i=0; i<g->n; i++){
		for(e=g->vertices[i]->neighborhood; e; e=e->next){
			int w = e->endPoint->number;
			if(w > i){
				// check whether edge already exists in edges list
				int edge_index = edgeUnique(e, edges);
				// if edge not in edges list yet, then add it, and corresponding graphs index it appears in
				if(edge_index == -1){
					dvint_realize(edges_graphs);
					vint_append(edges_graphs->data[edges->size], index);
					vector_append(edges, e);
				}else{ // if edge already exist in edges list
					// check the edge's appears graphs contains current graph index
					// if not, add graph index here
					if(vint_contain(edges_graphs->data[edge_index], index) == -1){
						vint_append(edges_graphs->data[edge_index], index);
					}
				}
			}
		}
	}
}

void generateEdges(Vector* graphDB, Vector* edges, DVectorInt* graphsp){
	struct Graph* g = NULL;
	for(int i=0; i<graphDB->size; i++){
		g = graphDB->data[i];
		scanEdgesSingleGraph(g, edges, graphsp, i);
	}
}

/* Print edges for debugging
 **/
void printEdges(Vector* edges){
	struct VertexList* e;
	for(int i=0; i<edges->size; i++){
		e = edges->data[i];
		printf("start: %s, end: %s, edge: %s\n", e->startPoint->label, e->endPoint->label, e->label);
	}
}

/* Print vertices for debugging
 **/
void printVertices(Vector* vertices){
	struct Vertex* v;
	printf("Vertices: ");
	for(int i=0; i<vertices->size; i++){
		v = vertices->data[i];
		printf("%s ", v->label);
	}
	printf("\n");
}

/* check whether s1 and s2 are the same
 **/
//int sameLabel(char* s1, char* s2){
int sameLabel(struct ShallowGraph* s1, struct ShallowGraph* s2){
	//if(strcmp(s1, s2) == 0){
	if(compareCanonicalStrings(s1, s2) == 0){
		return 1;
	}
	return 0;
}

/* Check whether the start and end point of the edge are in our frequent single vertices
 **/
int edgeValid(struct VertexList* e, Vector* frequent_single_vertices){
	int start = 0, end = 0;
	for(int i=0; i<frequent_single_vertices->size; i++){
		struct Vertex* v = frequent_single_vertices->data[i];
		if(strcmp(e->startPoint->label, v->label) == 0){
			start = 1;
		}
		if(strcmp(e->endPoint->label, v->label) == 0){
			end = 1;
		}
		if(start == 1 && end == 1){
			return 1;
		}
	}
	return 0;
}

/* Duplication check by using canonical string comparison
 **/
//int isDuplicateString(char* candidate_cs, Vector* cs_table){
int isDuplicateString(struct ShallowGraph* candidate_cs, Vector* cs_table){
	for(int i=0; i<cs_table->size; i++){
		//char* sg = cs_table->data[i];
		struct ShallowGraph* sg = cs_table->data[i];
		//if(strcmp(candidate_cs, sg) == 0){
		if(compareCanonicalStrings(candidate_cs, sg) == 0){
			return 1; // duplicated
		}
	}
	return 0; // not duplicated
}

/* check whether cs is in vector of candidates object
 **/
//int newIsDuplicateString(char* cs, Vector* candidates){
//int newIsDuplicateString(char* cs, Vector* candidates){
int newIsDuplicateString(struct ShallowGraph* cs, Vector* candidates){
	struct CandidateObjects* object = NULL;
	for(int i=0; i<candidates->size; i++){
		object = candidates->data[i];
		if(sameLabel(object->canonicalString, cs)){
			return i;
		}
	}
	return -1;
}

/* get the canonical string
 **/
char* getCanonicalString(struct Graph* tree, struct ShallowGraphPool* sgp){
	struct ShallowGraph* cs = canonicalStringOfTree(tree, sgp);
	char* canonical= canonicalStringToChar(cs);
	dumpShallowGraph(sgp, cs);
	return canonical;
}

/* check whether value v in vector vint
 **/
int vintIn(VectorInt* vint, int v){
	for(int i=0; i<vint->size; i++){
		if(vint->data[i] == v){
			return 1;
		}
	}
	return 0;
}

/* compare two ordered array, and output their intersection in O(n)
 **/
void orderedIntersection(VectorInt* v1, VectorInt* v2, VectorInt* intersection, int* i, int* j){
	for(; i[0]<v1->size && j[0]<v2->size;){
		if(v1->data[i[0]] == v2->data[j[0]]){
			vint_append(intersection, v1->data[i[0]]);
			++i[0];
			++j[0];
			orderedIntersection(v1, v2, intersection, i, j);
		}else if(v1->data[i[0]] < v2->data[j[0]]){
			++i[0];
			orderedIntersection(v2, v1, intersection, j, i);
		}else{
			++j[0];
			orderedIntersection(v1, v2, intersection, i, j);
		}
	}
}

/* get the intersection between two integer vector in O(n)
 **/
VectorInt* intersection(VectorInt* v1, VectorInt* v2){
	VectorInt* intersect = malloc(sizeof(VectorInt));
	vint_init(intersect);
	int i=0, j=0;
	orderedIntersection(v1, v2, intersect, &i, &j);
	/*
	for(int i=0; i<v1->size; i++){
		if(vintIn(v2, v1->data[i])){
			vint_append(intersect, v1->data[i]);
		}
	}*/
	return intersect;
}

/* get the intersection between the support vectors in sgs, by using index from parents
 **/
VectorInt* calculateIntersection(VectorInt* parents, Vector* sgs){
	VectorInt *intersect = NULL, *intersect_next = NULL;
	// get the intersection of first two support
	struct CandidateObjects* obj1 = sgs->data[parents->data[0]];
	struct CandidateObjects* obj2 = sgs->data[parents->data[1]];
	intersect = intersection(obj1->graphs, obj2->graphs);
	// iteratively calculate intersection
	for(int i=2; i<parents->size; i++){
		struct CandidateObjects* obj1 = sgs->data[parents->data[i]];
		intersect_next = intersection(intersect, obj1->graphs);
		vint_free(intersect);
		free(intersect);
		intersect = intersect_next;
	}
	// return the intersection of all support
	return intersect;
}

int binarySearch(char* cs, Vector* sgs){
	int first = 0, last = sgs->size-1;
	while(first <= last){
		int middle = (last+first) / 2;

		char* object = sgs->data[middle];

		if(strcmp(cs, object) == 0){
			///printf("found at location %d.\n", object->number);
			return middle;
		}else if(strcmp(cs, object) > 0){
			first = middle + 1;
		}else{
			last = middle - 1;
		}
	}
	//printf("Not found! %s is not present in the list.\n", cs);
	return -1;
}

//int binaryChecking(char* cs, Vector* sgs){
int binaryChecking(struct ShallowGraph* cs, Vector* sgs){
	int first = 0, last = sgs->size-1;
	while(first <= last){
		int middle = (last+first) / 2;

		struct CandidateObjects* object = sgs->data[middle];

		//if(strcmp(cs, object->canonicalString) == 0){  compareCanonicalStrings
		if(compareCanonicalStrings(cs, object->canonicalString) == 0){
			///printf("found at location %d.\n", object->number);
			return middle;
		//}else if(strcmp(cs, object->canonicalString) > 0){
		}else if(compareCanonicalStrings(cs, object->canonicalString) > 0){
			first = middle + 1;
		}else{
			last = middle - 1;
		}
	}
	//printf("Not found! %s is not present in the list.\n", cs);
	return -1;
}

/* Check whether each canonical string in sg are all in canonical string array sgs
 **/
VectorInt* existIn(Vector* sg, Vector* sgs){
	VectorInt* vint = malloc(sizeof(VectorInt));
	vint_init(vint);
	// for all the possible sub canonical string of tree canonical string
	for(int i=0; i<sg->size; i++){
		// check whether it is duplicated with any previous frequent graph
		//int result = newIsDuplicateString(sg->data[i], sgs);
		int result = binaryChecking(sg->data[i], sgs);
		// if it is duplicated with any previous frequent graph, get graph number
		if(result != -1){
			vint_append(vint, result);
		}else{ // if any sub canonical string is not in previous frequent graph
			vint_free(vint);
			free(vint);
			return NULL;
		}
	}
	// if all of the sub trees exist in previous level, calculate their support intersection
	//int size = 0;
	// if only one subtree ()
	/*if(vint->size == 1){
		obj = sgs->data[vint->data[0]];
		int size = obj->graphs->size;
		vint_free(vint);
		free(vint);
		return size;
	}else{*/
	VectorInt* intersect = calculateIntersection(vint, sgs);
	vint_free(vint);
	free(vint);
	return intersect;
	//}
}

/* Generate all the possible one vertex less subtree canonical string, by using graphs vertex removing
 * the vertex to remove must be a leaf, and also much be extendible from frequent single vertex and edges
 **/
void reducedShallowGraphsCanonicalString(Vector* reduced_graphs, struct Graph* g, struct GraphPool* gp, struct ShallowGraphPool* sgp, struct VertexPool* vp){
	struct Graph* g_reduced;
	for(int i=0; i<g->n; i++){
		if(isLeaf(g->vertices[i])){ // leaf and removable, all the leafs should be removable, as we extend from single vertex  && isRemovable(g->vertices[i], edges)
			g_reduced = createSGraph(gp, vp, g, i);
			struct ShallowGraph* sg = canonicalStringOfTree(g_reduced, sgp);
			//printShallowGraph(sg);
			//vector_append(reduced_graphs, canonicalStringToChar(sg));
			vector_append(reduced_graphs, sg);
			//dumpShallowGraph(sgp, sg);
			dumpGraph(gp, g_reduced);
		}
	}
}

/* Generate all the possible one vertex less subtree canonical string, from tree canonical string
 **/
void reducedShallowGraphsCS(Vector* reduced_sgs, char* cs){
	int left_ = -1, right_ = -1;
	VectorInt pairs;
	vint_init(&pairs);
	// get the range of leafs in canonical string into pairs
	int i = 0;
	for(; cs[i] != '\0'; i++){
		if(cs[i] == '('){
			left_ = i;
		}
		if(cs[i] == ')'){
			right_ = i;
			if(left_ != -1){
				vint_append(&pairs, left_);
				vint_append(&pairs, right_);
			}
			left_ = -1;
		}
	}
	// get all the sub canonical strings
	for(int j=0; j<pairs.size; j+=2){
		char* copy = (char*)malloc(i * sizeof(char));
		int size = 0;
		for(int k=0; cs[k] != '\0'; k++){
			// here use pairs.data[j] -1 to skip one whitespace
			if(k < pairs.data[j] - 1 || k > pairs.data[j+1]){
				copy[size++] = cs[k];
			}
		}
		copy[size] = '\0';
		printf("Reduced: %s\n", copy);
		vector_append(reduced_sgs, copy);
	}
}

/* free generated reduced canonical string
 **/
void dumpReducedGraphsCS(Vector* reduced, struct ShallowGraphPool* sgp){
		for(int j=0; j<reduced->size; j++){
			//free(reduced->data[j]);
			dumpShallowGraph(sgp, reduced->data[j]);
		}
		vector_free(reduced);
}

/* Create graph with single vertex
 **/
struct Graph* createSingleVertexGraph(struct GraphPool* gp, char* label){
	struct Graph* g = getGraph(gp);
	int v=0;
	setVertexNumber(g, 1);

	g->vertices[v] = getVertex(gp->vertexPool);
	g->vertices[v]->number = v;
	g->vertices[v]->id = v;
	g->vertices[v]->label = label;
	g->vertices[v]->neighborhood = NULL;

	return g;
}

/* Write frequent graph information to file pointed by fp
 **/
void writeGraph(int number, struct Graph* g, int frequency, FILE* fp){
	g->number = number;
	g->activity = frequency;
	printGraphAidsFormat(g, fp);
}

/* Write corresponding frequent graphs p for graphs g
 **/
void writeGraphs(DVectorInt* graphsps, FILE* fp2){
	for(int i=0; i<graphsps->size; i++){
		for(int j=0; j<graphsps->data[i]->size; j++){
			fprintf(fp2, "%i ", graphsps->data[i]->data[j]);
		}
		fprintf(fp2, "\n");
	}
}

/* generate frequent vertices from graph database
 **/
Vector* generateFrequentVertices(Vector* graphDB, int t, struct GraphPool* gp, DVectorInt* graphsps, FILE* fp){
	Vector vertices;
	vector_init(&vertices);
	Vector* frequent_vertices = malloc(sizeof(Vector));
	vector_init(frequent_vertices);
	DVectorInt graphsp;
	dvint_init(&graphsp);
	int number_of_frequent_graphs = -1;
	// generate all the vertices with different label
	generateVertices(graphDB, &vertices, &graphsp);
	// remove not frequent vertices
	for(int i=0; i<vertices.size; i++){
		if(graphsp.data[i]->size >= t){
			++number_of_frequent_graphs;
			for(int j=0; j<graphsp.data[i]->size; j++){
				vint_append(graphsps->data[j], number_of_frequent_graphs);
			}
			struct Vertex* v = vertices.data[i];
			struct Graph* g = createSingleVertexGraph(gp, v->label);
			writeGraph(number_of_frequent_graphs, g, graphsp.data[i]->size, fp);
			dumpGraph(gp, g);
			vector_append(frequent_vertices, vertices.data[i]);
		}
	}
	//garbage collection
	dvint_free(&graphsp);
	vector_free(&vertices);
	return frequent_vertices;
}

/* generate frequent edges from graph database
 **/
Vector* generateFrequentEdges(Vector* graphDB, Vector* edges_graphs, int t){
	Vector edges;
	vector_init(&edges);
	Vector* frequent_edges = malloc(sizeof(Vector));
	vector_init(frequent_edges);
	DVectorInt graphsp;
	dvint_init(&graphsp);
	generateEdges(graphDB, &edges, &graphsp);
	// remove not frequent edges
	for(int i=0; i<edges.size; i++){
		if(graphsp.data[i]->size >= t){
			vector_append(frequent_edges, edges.data[i]);
			vector_append(edges_graphs, vint_copy(graphsp.data[i]));
		}
	}
	// garbage collection
	dvint_free(&graphsp);
	vector_free(&edges);
	return frequent_edges;
}

/* free a vector of VectorInt
 **/
void freeVectorInt(Vector* vector_int){
	for(int i=0; i<vector_int->size; i++){
		vint_free(vector_int->data[i]);
		free(vector_int->data[i]);
	}
	vector_free(vector_int);
}

/* get the order of a candidate object by canonical string order
 **/
int orderedCanonicalStringVector(Vector *sgs, struct CandidateObjects *csg){
	struct CandidateObjects *csg_ = NULL;
	for(int i=0; i<sgs->size; i++){
		csg_ = sgs->data[i];
		//if(strcmp(csg->canonicalString, csg_->canonicalString) < 0){ // order canonical strings
			//return i;
			/*
			csg->number = i;
			for(int j=i+1; j<sgs->size; j++){
				csg_ = sgs->data[j];
				csg_->number = j;
			}*/
			//break;
		//}
		if(compareCanonicalStrings(csg->canonicalString, csg_->canonicalString) < 0){
			return i;
			break;
		}

	}
	//printf("index number: %i\n", csg->number);
	return csg->number;
}

/* transform a single edge to graph
 **/
struct Graph* edgeToGraph(struct VertexList* edge, struct GraphPool* gp){
	struct Graph* g = getGraph(gp);
	setVertexNumber(g, 2);
	// get vertices of edge
	for(int i=0; i<2; i++){
		g->vertices[i] = getVertex(gp->vertexPool);
		g->vertices[i]->number = i;
		g->vertices[i]->id = i;
	}
	g->vertices[0]->label = edge->startPoint->label;
	g->vertices[1]->label = edge->endPoint->label;
	// connect edge
	addEdgeBetweenVertices(0, 1, edge->label, g, gp);
	g->n = 2;
	g->m = 1;
	return g;
}

void printfOrder(VectorInt* graphs_order){
	printf("Order: ");
	for(int i=0; i<graphs_order->size; i++){
		printf("%i ", graphs_order->data[i]);
	}
	printf("\n");
}

/* generate CandidateObjects from a vector of edges
 **/
void generateGraphsFromEdges(Vector* frequent_edges, Vector* graphs, Vector* edges_appear_graphs, DVectorInt* graphsps, int number_of_frequent_graphs, struct GraphPool* gp, struct ShallowGraphPool* sgp, FILE* fp){
	struct VertexList* edge = NULL;
	struct CandidateObjects *object = NULL;
	for(int i=0; i<frequent_edges->size; i++){
		object = malloc(sizeof(struct CandidateObjects));
		++number_of_frequent_graphs;
		edge = frequent_edges->data[i];
		struct Graph* g = edgeToGraph(edge, gp);
		VectorInt* gs = edges_appear_graphs->data[i];
		// create objects
		object->characteristics = NULL;
		//object->canonicalString = getCanonicalString(g, sgp);
		object->canonicalString = canonicalStringOfTree(g, sgp);
		object->number = graphs->size;
		object->h = g;
		object->graphs = gs;
		int string_order = orderedCanonicalStringVector(graphs, object);
		//vint_insert(graphs_order, string_order, object->number);
		//printfOrder(graphs_order);
		//vector_append(graphs, object);
		vector_insert(graphs, string_order, object);
		// write to file
		for(int j=0; j<gs->size; j++){
			vint_append(graphsps->data[j], number_of_frequent_graphs);
		}
		writeGraph(number_of_frequent_graphs, g, gs->size, fp);
	}
}

int newBackCheckingSingle(Vector* pre_graphs, struct CandidateObjects* candidate, int t, struct GraphPool* gp, struct ShallowGraphPool* sgp, struct VertexPool* vp){
	Vector reduced;
	vector_init(&reduced);
	// remove one leaf from tree, and rearrange the number of the reduced tree
	reducedShallowGraphsCanonicalString(&reduced, candidate->h, gp, sgp, vp);
	// if all the reduced graphs are in previous level frequent graphs
	VectorInt* exist = existIn(&reduced, pre_graphs);
	dumpReducedGraphsCS(&reduced, sgp);
	if(exist != NULL){
		if(exist->size >= t){
			candidate->graphs = exist;
			return 1;
		}else{
			vint_free(exist);
			free(exist);
			return 0;
		}
	}
	return 0;
}

/* free the generated CandidateObjects in refine operator
 **/
void freeROgarbage(struct CandidateObjects* object, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	if(object->graphs){
		vint_free(object->graphs);
	}
	free(object->graphs);
	if(object->canonicalString){
		//free(object->canonicalString);
		dumpShallowGraph(sgp, object->canonicalString);
	}
	dumpGraph(gp, object->h);
	free(object);
}

int largerOrder(char* s1, char* s2){
	if(strcmp(s1, s2) <= 0){
		return 1;
	}
	return 0;
}

struct Vertex* extendibleVertex(struct VertexList* edge, struct Vertex* start){
	//if(sameLabel(edge->startPoint->label, start->label)){// && largerOrder(start->label, edge->endPoint->label)
	if(strcmp(edge->startPoint->label, start->label) == 0){
		return edge->endPoint;
	//}else if(sameLabel(edge->endPoint->label, start->label)){// && largerOrder(start->label, edge->startPoint->label)
	}else if(strcmp(edge->endPoint->label, start->label) == 0){// && largerOrder(start->label, edge->startPoint->label)
		return edge->startPoint;
	}else{
		return NULL;
	}
}

void duplicateCandidateElimination(struct CandidateObjects* object, struct CandidateObjects* pre_object, struct Vertex* candidateTree, Vector* candidates, char search, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	object->canonicalString = canonicalStringOfTree(object->h, sgp);
	if((search == 'p' && containsString(candidateTree, object->canonicalString)) || (search == 'y' && binaryChecking(object->canonicalString, candidates) != -1)){
		freeROgarbage(object, gp, sgp);
	}else{
		struct ShallowGraph* cs_copy = cloneShallowGraph(object->canonicalString, sgp);
		addToSearchTree(candidateTree, cs_copy, gp, sgp);
		object->parent = pre_object;
		int order = orderedCanonicalStringVector(candidates, object);
		vector_insert(candidates, order, object);
	}
}

/* Generate a vector of candidate objects, which contain the following formation
 * 1. graph h: struct Graph*, free by dumpGraph(gp, h)
 * 2. canonical string cs: char*, free by free(cs)
 * 3. which graphs it may appear: VectorInt*, free by vint_free(vint), free(vint)
 **/
int newRefineOperator(Vector* Graphs, Vector* candidates, Vector* frequent_edges, struct GraphPool* gp, struct ShallowGraphPool* sgp, char bd, char search_method, int t){
	// iterate all the graphs in S
	int count = 0;
	//struct Vertex* searchTree = searchTreeOneLevel(Graphs, gp, sgp);
	//struct ShallowGraph* strings = NULL;
	//struct Vertex* candidateTree = buildSearchTree(strings, gp, sgp); // create an empty prefix search tree
	struct Vertex* candidateTree = getVertex(gp->vertexPool);
	candidateTree->neighborhood = NULL;
	//candidateTree->neighborhood = NULL;
	for(int i=0; i<Graphs->size; i++){
		struct CandidateObjects* pre_object = Graphs->data[i];
		// extend over all the vertice of graph
		for(int v=0; v<pre_object->h->n; v++){
			struct Vertex* gver = pre_object->h->vertices[v];
			// for all the frequent edges
			for(int j=0; j<frequent_edges->size; j++){
				struct VertexList* e = frequent_edges->data[j];
				struct Vertex* ver = extendibleVertex(e, gver);
				// if edge is extendible with these two vertices
				if(ver){
					++count;
					// create a new candidate object
					struct CandidateObjects* object = malloc(sizeof(struct CandidateObjects));
					// get a copy of p from S
					object->h = cloneGraph(pre_object->h, gp);
					// initialize characteristics and graphs appear of candidate object
					object->characteristics = NULL;
					object->graphs = NULL;
					object->canonicalString = NULL;
					object->number = candidates->size;
					// extend p by adding a single vertex with the lowest vertex number
					addSingleVertex(object->h, object->h->vertices[v], ver, gp, e->label);
					// check whether all subtrees of candidate exist in previous level if using level wise method
					if((bd == 'b' && newBackCheckingSingle(Graphs, object, t, gp, sgp, gp->vertexPool)) || bd == 'd'){
						duplicateCandidateElimination(object, pre_object, candidateTree, candidates, search_method, gp, sgp);
						//object->canonicalString = getCanonicalString(object->h, sgp);
						/*
						object->canonicalString = canonicalStringOfTree(object->h, sgp);
						//if(binaryChecking(object->canonicalString, candidates) != -1){
						if(containsString(candidateTree, object->canonicalString)){
							freeROgarbage(object, gp, sgp);
						}else{
							struct ShallowGraph* cs_copy = cloneShallowGraph(object->canonicalString, sgp);
							addToSearchTree(candidateTree, cs_copy, gp, sgp);
							object->parent = pre_object;
							//vector_append(candidates, object);
							int order = orderedCanonicalStringVector(candidates, object);
							vector_insert(candidates, order, object);
						}*/
					}else{
						freeROgarbage(object, gp, sgp);
					}
				}
			}
		}
	}
	dumpSearchTree(gp, candidateTree);
	return count;
}

/* free a single CandidateObjects
 **/
void freeCandidateObject(struct CandidateObjects *object, struct GraphPool* gp, struct ShallowGraphPool* sgp, char bd){
	dumpGraph(gp, object->h);
	if(bd == 'b'){ // we keep all the canonical string for depth first search
		//free(object->canonicalString);
		dumpShallowGraph(sgp, object->canonicalString);
		if(object->graphs != NULL){
			vint_free(object->graphs);
			free(object->graphs);
		}
	}
	if(object->characteristics){
		int i=0;
		while(i < object->characteristics->size){
			VectorCharacteristics* vc = object->characteristics->data[i];
			vc_free(vc);
			free(vc);
			++i;
		}
		vector_free(object->characteristics);
		free(object->characteristics);
	}
	/*
	if(object->graphs != NULL){
		vint_free(object->graphs);
		free(object->graphs);
	}*/
	free(object);
}

/* Generate frequent 3rd order graphs, each graph contains the following information
 * 1. struct Graph*, frequent 3rd order graph
 * 2. Vector* characteristics, a vector of VectorCharacteristics objects
 * 3. VectorInt*, a integer vector contains information of which graphs the candidate appears by graph index
 **/
void newGenerateThirdOrderGraphs(Vector* graphs_of_edges, DVectorInt* graphsp, Vector* graphs, Vector* frequent_edges, Vector* graphDB, Vector* objects, struct GraphPool* gp, struct ShallowGraphPool* sgp, char bd, char search_method, int t, char iterative){
	Vector candidates;
	vector_init(&candidates);
	int count = newRefineOperator(graphs_of_edges, &candidates, frequent_edges, gp, sgp, bd, search_method, t);
	printf("Number of extended candidates: %i\n", count);
	printf("Number of possible frequent candidates: %i\n", candidates.size);

	for(int i=0; i<candidates.size; i++){
		struct CandidateObjects* object = candidates.data[i];
		int frequency = 0;
		dvint_realize(graphsp);
		object->characteristics = NULL;
		if(iterative == 'f'){
			object->characteristics = malloc(sizeof(Vector));
			vector_init(object->characteristics);
		}
		// for all the graph candidate may appears
		for(int j=0; j<object->graphs->size; j++){ // for all the graphs, candidate's may appears in
			//VectorCharacteristics* vc = NULL;
			int graph_index = object->graphs->data[j];
			struct Graph* g = graphDB->data[graph_index];
			if(iterative == 'f'){
				VectorCharacteristics* vc = malloc(sizeof(VectorCharacteristics));
				vc_init(vc);
				if(subtreeIsomorphic(vc, object->h, g, objects->data[graph_index])){
					vector_append(object->characteristics, vc);
					vint_append(graphsp->data[graphs->size], graph_index);
					++frequency;
				}else{
					vc_free(vc);
					free(vc);
				}
			}else if(subGraphIsoCheck(object->h, g, objects->data[graph_index], gp, sgp)){
				vint_append(graphsp->data[graphs->size], graph_index);
				++frequency;
			}
		}
		if(frequency >= t){
			vint_free(object->graphs);
			free(object->graphs);
			object->graphs = graphsp->data[graphs->size];
			vector_append(graphs, object);
		}else{
			freeCandidateObject(object, gp, sgp, 'b');
			dvint_clean(graphsp);
		}
	}
	vector_free(&candidates);
	printf("Number of frequent trees with order 3: %i\n", graphs->size);
}

/* free a vector of CandidateObjects
 **/
void freeVectorGraphs(struct GraphPool* gp, struct ShallowGraphPool* sgp, Vector* graphs){
	for(int i=0; i<graphs->size; i++){
		freeCandidateObject(graphs->data[i], gp, sgp, 'b');
	}
	vector_free(graphs);
}

/* free generated frequent vertices and edges
 **/
void freeVerticesEdges(Vector* frequent_vertices, Vector* frequent_edges){
	vector_free(frequent_vertices);
	free(frequent_vertices);
	vector_free(frequent_edges);
	free(frequent_edges);
}

/* process file writing objects at the end of execution
 **/
void writeFile(DVectorInt* graphsps, FILE** file){
	// write corresponding frequent graphs p for graphs g
	writeGraphs(graphsps, file[1]);
	// add end symbol to file and close file object
	fprintf(file[0], "$");
	fclose(file[0]);
	fclose(file[1]);
	// free graph writing objects
	dvint_free(graphsps);
	free(file);
}

/* initialize the file writing information
 **/
FILE** fileWriting(Vector* GraphDB, DVectorInt* graphsps){
	FILE **file = malloc(2*sizeof(FILE*));
	for(int i=0; i<GraphDB->size; i++){
		dvint_realize(graphsps);
		vint_append(graphsps->data[i], ((struct Graph*)GraphDB->data[i])->number);
		vint_append(graphsps->data[i], ((struct Graph*)GraphDB->data[i])->activity);
	}
	file[0] = fopen("./frequent_graphs.txt", "w+");
	file[1] = fopen("./graphs.txt", "w+");
	return file;
}

/* get the sequence number in candidate's parent's support
 **/
int getPAppGraphCha(int graph_index, struct CandidateObjects* object){
	struct CandidateObjects* parent = object->parent;
	int i=0;
	for(; i<parent->graphs->size;){
		if(parent->graphs->data[i] == graph_index){
			break;
		}
		++i;
	}
	return i;
}

// free one level graphs, graph numbers and characteristics
void freeOneLevelObjects(Vector* graphs, DVectorInt* graphsp, DVectorInt* graphsps, struct GraphPool* gp, struct ShallowGraphPool* sgp, int* number_of_frequent_graphs, FILE* fp){
	for(int i=0; i<graphs->size; i++){
		struct CandidateObjects *obj = graphs->data[i];
		++(*number_of_frequent_graphs);
		for(int j=0; j<graphsp->data[i]->size; j++){
			vint_append(graphsps->data[j], *number_of_frequent_graphs);
		}
		// write frequent graphs in S{l} to file
		writeGraph(*number_of_frequent_graphs, obj->h, obj->graphs->size, fp);
		freeCandidateObject(graphs->data[i], gp, sgp, 'b');
	}
	free(graphsp->data);
}

/* filter out non frequent candidates for both level wise and frequent pattern method
 **/
void filterNonFrequentCandidates(Vector* GraphDB, Vector* objects, Vector* all_graphs_cs, struct Vertex* search_tree, Vector* candidates, DVectorInt* graphsp_next, Vector* graphs_next, char flag, char db, int t, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	// iterate all the candidates
	for(int j=0; j<candidates->size; j++){
		struct CandidateObjects* object = candidates->data[j];
		struct Graph* candidate = object->h;
		int frequency = 0;
		dvint_realize(graphsp_next);
		object->characteristics = NULL;
		if(flag == 'f'){ // if use iterative isomorphic
			object->characteristics = malloc(sizeof(Vector));
			vector_init(object->characteristics);
		}
		VectorInt* appears_in_graphs = NULL;
		if(db == 'b'){ // if use level wise method, used intersected support
			appears_in_graphs = object->graphs;
		}else{	// if use frequent pattern method, use parent'd support
			appears_in_graphs = object->parent->graphs;
		}
		// for all the graphs, candidate may appear
		for(int index=0; index<appears_in_graphs->size; index++){
			int graph_index = appears_in_graphs->data[index];
			struct Graph* g = GraphDB->data[graph_index];
			if(flag == 'f'){
				// b is always the last vertex number of candidate, and we extend b from vertex a
				struct Vertex* b = candidate->vertices[candidate->n - 1];
				struct Vertex* a = b->neighborhood->endPoint;
				VectorCharacteristics* vc = malloc(sizeof(VectorCharacteristics));
				vc_init(vc);
				// CC_pre will keep the characteristics from p w.r.t each graph from GraphDB
				int cha = -1;
				if(db == 'b'){ // if use level wise, intersect support is not same as parent's support
					cha = getPAppGraphCha(graph_index, object);
				}else{	// if use frequent pattern, then it is just the index
					cha = index;
				}
				if(fastSubGraphIsoCheck(candidate, g, objects->data[graph_index], object->parent->characteristics->data[cha], vc, a->number, b->number)){
					vector_append(object->characteristics, vc);
					vint_append(graphsp_next->data[graphsp_next->size-1], graph_index);
					++frequency;
				}else{
					vc_free(vc);
					free(vc);
				}
			}else if(subGraphIsoCheck(candidate, g, objects->data[graph_index], gp, sgp)){
				vint_append(graphsp_next->data[graphsp_next->size-1], graph_index);
				++frequency;
			}
		}
		// if candidate appears in more than t graphs in GraphDB
		if(frequency >= t){
			if(db == 'b'){
				vint_free(object->graphs);
				free(object->graphs);
			}else{
				vector_append(all_graphs_cs, object->canonicalString);
				struct ShallowGraph* cstring = cloneShallowGraph(object->canonicalString, sgp);
				addToSearchTree(search_tree, cstring, gp, sgp);
			}
			object->graphs = graphsp_next->data[graphsp_next->size-1];
			// add from candidate set C_next to S_next
			vector_append(graphs_next, object);
		}else{
			freeCandidateObject(object, gp, sgp, 'b'); // free the canonical string, so use 'b' here
			// free graph numbers
			dvint_clean(graphsp_next);
		}
	}
}

int generateVerticesEdges(Vector *GraphDB, int t, struct GraphPool* gp, struct ShallowGraphPool* sgp, DVectorInt* graphsps, FILE* fp, Vector** frequent_vertices, Vector** frequent_edges, Vector* graphs_of_edges){
	frequent_vertices[0] = generateFrequentVertices(GraphDB, t, gp, graphsps, fp);
	printf("Number of frequent vertices: %i\n", frequent_vertices[0]->size);
	int number_of_frequent_graphs = frequent_vertices[0]->size-1;
	Vector edges_appears_graphs;
	vector_init(&edges_appears_graphs);
	frequent_edges[0] = generateFrequentEdges(GraphDB, &edges_appears_graphs, t);
	printf("Number of frequent edges: %i\n", frequent_edges[0]->size);
	generateGraphsFromEdges(frequent_edges[0], graphs_of_edges, &edges_appears_graphs, graphsps, number_of_frequent_graphs, gp, sgp, fp);
	number_of_frequent_graphs += frequent_edges[0]->size;
	vector_free(&edges_appears_graphs);
	return number_of_frequent_graphs;
}

int generateThirdOrder(int order, Vector* graphs_of_edges, DVectorInt* graphsp, Vector* graphs, Vector* frequent_vertices, Vector* frequent_edges, Vector* GraphDB, Vector* objects, struct GraphPool* gp, struct ShallowGraphPool* sgp, char db, char search_method, int t, char flag, DVectorInt* graphsps, FILE** file){
	if(order >= 3){
		vector_init(graphs);
		newGenerateThirdOrderGraphs(graphs_of_edges, graphsp, graphs, frequent_edges, GraphDB, objects, gp, sgp, db, search_method, t, flag);
		return 1;
	}else{
		freeVerticesEdges(frequent_vertices, frequent_edges);
		freeVectorGraphs(gp, sgp, graphs_of_edges);
		writeFile(graphsps, file);
		free(graphsp->data);
		return 0;
	}
}

/* Frequent subgraph mining from graph database GraphDB, with level wise method
 **/
void apFrequentSubGraphMining(Vector* GraphDB, Vector* objects, struct GraphPool* gp, struct ShallowGraphPool* sgp, int t, int order, char flag, char search_method){
	// keep frequent graphs p for graph g
	DVectorInt graphsps;
	dvint_init(&graphsps);
	FILE** file = fileWriting(GraphDB, &graphsps);
	// get frequent vertices
	Vector *frequent_vertices, *frequent_edges;
	Vector graphs_of_edges;
	vector_init(&graphs_of_edges);
	int number_of_frequent_graphs = generateVerticesEdges(GraphDB, t, gp, sgp, &graphsps, file[0], &frequent_vertices, &frequent_edges, &graphs_of_edges);
	// generate 3rd order frequent graphs
	Vector graphs;
	// hold number of corresponding graphs g for frequent graph p
	DVectorInt graphsp;
	dvint_init(&graphsp);
	if(!generateThirdOrder(order, &graphs_of_edges, &graphsp, &graphs, frequent_vertices, frequent_edges, GraphDB, objects, gp, sgp, 'b', search_method, t, flag, &graphsps, file)){
		return;
	}
	int size = 3;
	freeVectorGraphs(gp, sgp, &graphs_of_edges);
	// while we still get frequent graph in S
	while(graphs.size && size < order){
		// keep the possible graphs for next level
		Vector graphs_next;
		vector_init(&graphs_next);
		// candidates hold all the possible non duplicate candidates from graph p in S
		Vector candidates;
		vector_init(&candidates);
		clock_t tic = clock();
		int count = newRefineOperator(&graphs, &candidates, frequent_edges, gp, sgp, 'b', search_method, t);
		clock_t toc = clock();
		printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		printf("Number of extended candidates: %i\n", count);
		printf("Number of possible frequent candidates: %i\n", candidates.size);
		// hold number of corresponding graphs g for frequent graph p
		DVectorInt graphsp_next;
		dvint_init(&graphsp_next);
		filterNonFrequentCandidates(GraphDB, objects, NULL, NULL, &candidates, &graphsp_next, &graphs_next, flag, 'b', t, gp, sgp);
		// free previous level graphs, and graph canonical string
		freeOneLevelObjects(&graphs, &graphsp, &graphsps, gp, sgp, &number_of_frequent_graphs, file[0]);
		graphsp = graphsp_next;
		// free last level graphs pointers, candidates, flags
		vector_free(&graphs);
		vector_free(&candidates);
		// redirect S to S_next
		graphs = graphs_next;
		printf("Number of frequent trees with order %i: %i\n", ++size, graphs.size);
	}
	// free last level graphs, graph numbers and characteristics
	freeOneLevelObjects(&graphs, &graphsp, &graphsps, gp, sgp, &number_of_frequent_graphs, file[0]);
	writeFile(&graphsps, file);
	vector_free(&graphs);
	// free frequent vertices and edges
	freeVerticesEdges(frequent_vertices, frequent_edges);
}

/* write graph for frequent pattern
 **/
void writeGraphss(int* number_of_graphs, Vector* graphs, FILE* fp){
	struct CandidateObjects* obj = NULL;
	for(int i=0; i<graphs->size; i++){
		obj = graphs->data[i];
		writeGraph(*number_of_graphs, obj->h, obj->graphs->size, fp);
		*number_of_graphs += 1;
	}
}

/* copy canonical string from vector of CanididateObjects
 **/
void copyCS(Vector* graphs_cs, Vector* graphs){
	for(int i=0; i<graphs->size; i++){
		struct CandidateObjects* obj = graphs->data[i];
		vector_append(graphs_cs, obj->canonicalString);
	}
}

/* Remove 1. non frequent pattern graphs.
 * 		  2. pattern graphs duplicate with already generated pattern graphs in previous recursive call
 **/
void removeNonFrequentCandidates(char flag, Vector* all_graphs_cs, struct Vertex* search_tree, Vector* candidates, Vector* graphs_next, DVectorInt* graphsp, float t, Vector* GraphDB, Vector* objects, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	// keep non duplicate candidates
	Vector non_duplicate;
	vector_init(&non_duplicate);
	// if candidate is duplicate with any already generated candidates, remove it
	for(int i=0; i<candidates->size; i++){
		struct CandidateObjects* obj = candidates->data[i];
		if(containsString(search_tree, obj->canonicalString)){
			freeCandidateObject(obj, gp, sgp, 'b');
		}else{
			vector_append(&non_duplicate, obj);
		}
		/*
		if(isDuplicateString(obj->canonicalString, all_graphs_cs)){
		//if(binarySearch(obj->canonicalString, all_graphs_cs) != -1){
			freeCandidateObject(obj, gp, sgp, 'b'); // also free canonical string
		}else{ // else append to non duplicated candidates
			vector_append(&non_duplicate, obj);
		}*/
	}
	filterNonFrequentCandidates(GraphDB, objects, all_graphs_cs, search_tree, &non_duplicate, graphsp, graphs_next, flag, 'd', t, gp, sgp);
	vector_free(&non_duplicate);
}

/* Once we don't have frequent graphs, then the iterative call will stop
 **/
void fpGraphGrowth(char flag, char search_method, int* number_of_graphs, Vector* all_graphs_cs, struct Vertex* search_tree, Vector* graphs, DVectorInt* graphsp, int t, int order, Vector* frequent_single_vertices, Vector* edges, Vector* GraphDB, Vector* objects, struct GraphPool* gp, struct ShallowGraphPool* sgp, FILE* fp){
	// for all frequent graphs p in S
	for(int i=0; i<graphs->size && ((struct CandidateObjects*)graphs->data[0])->h->n<order; i++){
		// candidates hold all the possible non duplicate candidates from graph p in S
		Vector candidates;
		vector_init(&candidates);
		// hold only one single frequent graph from graphs
		Vector single_graphs;
		vector_init(&single_graphs);
		vector_append(&single_graphs, graphs->data[i]);
		// get non duplicate candidates for a single p
		newRefineOperator(&single_graphs, &candidates, edges, gp, sgp, 'd', search_method, t);
		vector_free(&single_graphs);
		// remove non frequent candidates
		Vector graphs_next;
		vector_init(&graphs_next);
		removeNonFrequentCandidates(flag, all_graphs_cs, search_tree, &candidates, &graphs_next, graphsp, t, GraphDB, objects, gp, sgp);
		// write graphs in file
		writeGraphss(number_of_graphs, &graphs_next, fp);
		// recursively call fpGraphGrowth
		fpGraphGrowth(flag, search_method, number_of_graphs, all_graphs_cs, search_tree, &graphs_next, graphsp, t, order, frequent_single_vertices, edges, GraphDB, objects, gp, sgp, fp);
		for(int i=0; i<graphs_next.size; i++){
			freeCandidateObject(graphs_next.data[i], gp, sgp, 'd'); // not free canonical string
		}
		vector_free(&graphs_next);
		vector_free(&candidates);
	}
}

/* Get the order of a tree from it's canonical string representation
 **/
int canonicalStringOrder(char* cs){
	int order = 1;
	for(int i=0; cs[i]!='\0'; i++){
		if(cs[i] == '('){
			++order;
		}
	}
	free(cs);
	return order;
}

/* Output statics information of frequent pattern mining method, and also get writing information, then free it
 **/
void outputStatics(Vector* graphs_cs, DVectorInt* graphsp, DVectorInt* graphsps){
	VectorInt* statics = malloc(sizeof(VectorInt));
	vint_init(statics);
	for(int i=0; i<graphs_cs->size; i++){
		for(int j=0; j<graphsp->data[i]->size; j++){
			vint_append(graphsps->data[j], i);
		}
		int order = canonicalStringOrder(canonicalStringToChar(graphs_cs->data[i]));
		if(statics->size < order - 2){
			vint_append(statics, 1);
		}else{
			statics->data[order-3] += 1;
		}
		//free(graphs_cs->data[i]);
	}
	for(int i=1; i<statics->size; i++){
		printf("Number of frequent trees with order %i: %i\n", i+3, statics->data[i]);
	}
	vint_free(statics);
	free(statics);
}

/* build a search tree from a list of canonical strings without influence the original canonical string by copy
 **/
void buildsSearchTree(Vector* all_graphs_cs, struct Vertex* search_tree, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	for(int i=0; i<all_graphs_cs->size; i++){
		struct ShallowGraph* cstring = cloneShallowGraph(all_graphs_cs->data[i], sgp);
		addToSearchTree(search_tree, cstring, gp, sgp);
	}
}

/* Frequent subgraph mining from graph database GraphDB, with FP-growth
 * */
void fpFrequentSubGraphMining(Vector* GraphDB, Vector* objects, struct GraphPool* gp, struct ShallowGraphPool* sgp, int t, int order, char flag, char search_method){
	// keep frequent graphs p for graph g
	DVectorInt graphsps;
	dvint_init(&graphsps);
	FILE** file = fileWriting(GraphDB, &graphsps);
	// get frequent vertices, frequent edges, frequent graphs of edges
	Vector *frequent_vertices, *frequent_edges;
	Vector graphs_of_edges;
	vector_init(&graphs_of_edges);
	int number_of_frequent_graphs = generateVerticesEdges(GraphDB, t, gp, sgp, &graphsps, file[0], &frequent_vertices, &frequent_edges, &graphs_of_edges);
	// generate 3rd order frequent graphs
	Vector graphs;
	Vector all_graphs_cs;
	// hold number of corresponding graphs g for frequent graph p
	DVectorInt graphsp;
	dvint_init(&graphsp);
	if(!generateThirdOrder(order, &graphs_of_edges, &graphsp, &graphs, frequent_vertices, frequent_edges, GraphDB, objects, gp, sgp, 'b', search_method, t, flag, &graphsps, file)){
		return;
	}else{
		vector_init(&all_graphs_cs);
		copyCS(&all_graphs_cs, &graphs);
	}
	freeVectorGraphs(gp, sgp, &graphs_of_edges);
	// write graph into file
	writeGraphss(&number_of_frequent_graphs, &graphs, file[0]);
	// get prefix tree
	struct Vertex* search_tree = getVertex(gp->vertexPool);
	search_tree->neighborhood = NULL;
	buildsSearchTree(&all_graphs_cs, search_tree, gp, sgp);
	// frequent pattern growth method
	fpGraphGrowth(flag, search_method, &number_of_frequent_graphs, &all_graphs_cs, search_tree, &graphs, &graphsp, t, order, frequent_vertices, frequent_edges, GraphDB, objects, gp, sgp, file[0]);
	outputStatics(&all_graphs_cs, &graphsp, &graphsps);
	for(int i=0; i<graphs.size; i++){
		freeCandidateObject(graphs.data[i], gp, sgp, 'd');
	}
	dvint_free(&graphsp);
	// write corresponding frequent graphs p for graphs g
	writeFile(&graphsps, file);
	// free non used vectors
	vector_free(&all_graphs_cs);
	vector_free(&graphs);
	freeVerticesEdges(frequent_vertices, frequent_edges);
	dumpSearchTree(gp, search_tree);
}
