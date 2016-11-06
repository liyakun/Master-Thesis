#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include "./subgraphIsomorphism.h"
#include "./graphPrinting.h"
#include "./loading.h"
#include "./outerplanar.h"
#include "./subtreeIsomorphism.h"
#include "./listComponents.h"
#include "./listSpanningTrees.h"
#include "./frequentsubgraph.h"

/* check whether 'neighbor' is really w neighbor
 **/
int checkNeighbor(struct Vertex* w, struct Vertex* neighbor){
	for(struct VertexList* edge=w->neighborhood; edge!=NULL; edge=edge->next){
		// if the end points contain neighbor, return 1
		if(edge->endPoint->id == neighbor->id){
			return 1;
		}
	}
	// if the w does not have this neighbor
	return 0;
}

/* Print the edges of shallow graph for debugging
 **/
void printShallowTreeEdges(struct ShallowGraph* t){
	struct VertexList* e = NULL;
	for(e=t->edges; e; e=e->next){
		printf("(%i, %i) ", e->startPoint->id, e->endPoint->id);
	}
}

/* Print structure IsoCharacteristics for debugging
 **/
void printCharacteristic(struct IsoCharacteristics* characteristic){
	printf("(start point: %i", characteristic->start);
	printf(", tree edges:");
	printShallowTreeEdges(characteristic->t);
	printf(", end point: %i", characteristic->end);
	printf(", w: %i)\n", characteristic->w);
}

/* Initialize all attributes of structure IsoCharacteristics
 **/
void initCharacteristic(struct IsoCharacteristics* cha){
	cha->start = -1;
	cha->end = -1;
	cha->t = NULL;
	cha->w = -1;
}

int exist(struct ShallowGraph* t, int w, int c){
	for(struct VertexList* edge = t->edges; edge; edge=edge->next){
		if(((edge->startPoint->id == w) && (edge->endPoint->id == c))
				||((edge->startPoint->id == c) && (edge->endPoint->id == w))){
			return 1;
		}
	}
	return 0;
}

/* Get children of w, we have the post order of the vertex in a local spanning tree
 * 1. the children of a vertex is all the vertex post order left
 * 2. the children should also be the neighbor
 **/
void getChild(Vector* CT, struct Vertex* w, VectorInt* post_order, struct Graph* g, struct ShallowGraph* t){
	for(int i=0; post_order->data[i] != w->id; i++){
		//if(checkNeighbor(w, g->vertices[post_order->data[i]])){
		//	vector_append(CT, g->vertices[post_order->data[i]]);
		//}
		if(exist(t, w->id, post_order->data[i])){
			vector_append(CT, g->vertices[post_order->data[i]]);
		}
	}
}

/* get the neighbors of Vertex u from graph H
 **/
void getNeighbors(int u, int y, Vector* neighbors, struct Graph* H){
	// iterate all the neighborhood of vertex u
	for(struct VertexList* edge=H->vertices[u]->neighborhood; edge!=NULL; edge=edge->next){
		// if the id of y and u is different, we filter neighbor y of u, as required by the algorithm
		if(H->vertices[y]->id != H->vertices[u]->id){
			if(edge->endPoint->id != H->vertices[y]->id){ // filter y out from neighbors vertex array
				vector_append(neighbors, edge->endPoint);
			}
		}else{ // if y==u, keep all the neighbors
			vector_append(neighbors, edge->endPoint);
		}
	}
}

void getNeighbor(int u, int y, Vector* neighbors, struct Graph* H){
	for(struct VertexList* edge = H->vertices[u]->neighborhood; edge!=NULL; edge=edge->next){
		if(u!=y){
			if(edge->endPoint->number != y){
				vector_append(neighbors, edge->endPoint);
			}
		}else{
			vector_append(neighbors, edge->endPoint);
		}
	}
}

/* Check perfect bipartite matching for all the neighbors iteratively
 **/
int bipartiteMatrixMatching(int vu, int* seen, int* matchR, int size_theta, int (*bpGraph)[size_theta])
{
	// iterate all the possible c in CT&&CTPrime
	for (int v = 0; v < size_theta; v++){
		if (bpGraph[vu][v] && !seen[v]){
			seen[v] = 1;
			// if vertex has not been matched, or we are able to find a new match for it
			if (matchR[v] < 0 || bipartiteMatrixMatching(matchR[v], seen, matchR, size_theta, bpGraph)){
				matchR[v] = vu;
				return 1;
			}
		}
	}
	return 0;
}

/* Check the matching of neighbors of vertex u with c in CTheta in B
 **/
int checkMatching(Vector* b_vertices, int u, int y, struct Vertex* w, int num_neighbors, int size_theta, int (*bi_matrix)[size_theta], struct Graph* H){
	//for all the neighbors
	if(y != u){
		if(H->vertices[u]->neighborhood->endPoint->id == H->vertices[y]->id){
			//if the bipartite graph is empty and y is the only neighbor of u, and u have same label with w
			if(b_vertices->size==0 && H->vertices[u]->neighborhood->next == NULL && strcmp(H->vertices[u]->label, w->label) == 0){
				return 1;
			}
		}
	}
	if(num_neighbors==0 || size_theta==0){
		return 0;
	}
	int matchR[size_theta];
	memset(matchR, -1, sizeof(matchR));
	for(int vu=0; vu<num_neighbors; vu++){
		int seen[size_theta];
		memset(seen, 0, sizeof(seen));
		if (!bipartiteMatrixMatching(vu, seen, matchR, size_theta, bi_matrix))
			return 0; //once we can not find matching for any neighbor, return 0
	}
	return 1; //we find matching for all neighbors
}

/* Transfer an integer array to vector integer array, also transfer vertex number to vertex id
 **/
void intArrayToVint(VectorInt* t, int* post_order, struct Graph* g){
	for(int i=0; i<g->n; i++){
		vint_append(t, g->vertices[post_order[i]]->id);
	}
}

/* Get the post order of spanning tree t, rooted vertex v
 **/
VectorInt* getSpanningTreePostOrder(struct ShallowGraph* t, struct GraphPool* gp, struct Vertex* v){
	struct Graph* tmp_g = shallowGraphToGraph(t, gp);
	int* post_order = (int*)getPostorder(tmp_g, getNumById(tmp_g, v->id));
	VectorInt* it = malloc(sizeof(VectorInt));
	vint_init(it);
	intArrayToVint(it, post_order, tmp_g);
	free(post_order);
	dumpGraph(gp, tmp_g);
	return it;
}

/* Construct a IsoCharacteristics structure to tmp
 **/
void constructIsoCharacteristics(struct IsoCharacteristics* tmp, int start, int end, struct Vertex* c, struct ShallowGraph* t){
	tmp->start = start;
	tmp->end = end;
	tmp->w = c->id;
	tmp->t = t;
}

/* Get the vertices of the bipartite graph, we already remove y from neighbors here if y!=u
 **/
int get_c_u_primes(Vector* bi_c_u, struct ShallowGraph* t, struct ShallowGraph* t_prime, Vector* CT, Vector* CTPrime, int u, Vector* neighbors, VectorCharacteristics* CC){
	int num_neighbor = 0;
	for(int i=0; i<neighbors->size; i++){
		//check the c in CT
		int flag = 0;
		for(int m=0; m<CT->size; m++){
			struct IsoCharacteristics* tmp = (struct IsoCharacteristics*)malloc(sizeof(struct IsoCharacteristics));
			constructIsoCharacteristics(tmp, u, ((struct Vertex*)neighbors->data[i])->number, CT->data[m], t);
			if(vc_contain(CC, tmp)){
				vector_append(bi_c_u, neighbors->data[i]);
				vector_append(bi_c_u, CT->data[m]);
				flag = 1;
			}
			free(tmp);
		}
		//check the c in CTPrime
		for(int n=0; n<CTPrime->size; n++){
			struct IsoCharacteristics* tmp = (struct IsoCharacteristics*)malloc(sizeof(struct IsoCharacteristics));
			constructIsoCharacteristics(tmp, u, ((struct Vertex*)neighbors->data[i])->number, CTPrime->data[n], t_prime);
			if(vc_contain(CC, tmp)){
				vector_append(bi_c_u, neighbors->data[i]);
				vector_append(bi_c_u, CTPrime->data[n]);
				flag = 1;
			}
			free(tmp);
		}
		// if get a neighbor and also c from CC
		if(flag == 1){
		  num_neighbor++;
		}
	}
	// return the number of neighbors in CC
	return num_neighbor;
}

/* Join vertices in CT and CTPrime together to construct bipartite matrix
 **/
void joinCtCtprime(Vector* CTheta, Vector* CT, Vector* CTPrime){
	int i=0;
	while(i<CT->size){
		vector_append(CTheta, CT->data[i]);
		++i;
	}
	i=0;
	while(i<CTPrime->size){
		vector_append(CTheta, CTPrime->data[i]);
		++i;
	}
}

/* Get the index of the vertex in CTheta
 **/
int getIdx(Vector* CTheta, struct Vertex* c){
	int i=0;
	// here we should always find the index of c, as we get it from CT and CTPrime
	while(i< CTheta->size && ((struct Vertex*)CTheta->data[i])->id != c->id){
		++i;
	}
	return i;
}

/* Get the edge label between vertex neighbor and vertex v in graph g
 **/
char* getEdgeLabel(struct Vertex* neighbor, struct Vertex* v, struct Graph* g){
	struct VertexList* e;
	for(e=g->vertices[v->id]->neighborhood; e; e=e->next){
		if(e->endPoint->id == neighbor->id){
			return e->label;
		}
	}
	printf("ERROR: NOT FIND VERTEX PAIR\n");
	return NULL;
}

char* getEdgeLabelH(struct Vertex* neighbor, struct Vertex* u){
	struct VertexList* e;
	for(e=u->neighborhood; e; e=e->next){
		if(e->endPoint->number == neighbor->number){
			return e->label;
		}
	}
	printf("ERROR: NOT FIND VERTEX PAIR\n");
	return NULL;
}

/* Check whether the edge label in the pattern graph and the graph G is the same for specified vertex
 **/
int checkEdgeLabel(struct Vertex* neighbor, struct Vertex* c, int u, struct Vertex* w, struct Graph* H, struct Graph* G){
	if(strcmp(getEdgeLabelH(neighbor, H->vertices[u]), getEdgeLabel(c, w, G)) == 0){
		return 1;
	}
	return 0;
}

/* Construct a bipartite matrix for check a perfect matching between neighbors and c
 **/
void constructBiMatrix(int size_theta, int (*bi_matrix)[size_theta], Vector* bi_edges, Vector* CTheta, int u, struct Vertex* w, struct Graph* H, struct Graph* G){
	int j=0;
	// iterate all the vertex in bi_edges, which contains pairs (neighbor, c), increase 2 to index all neighbors
	for(int i=0; i < bi_edges->size; i+=2){
		// get the index of the c in position (i+1)
		int c_idx = getIdx(CTheta, (struct Vertex*)bi_edges->data[i+1]);
		// check the label of edges between vertex pairs (neighbor, u) and (c, w) in graph H and G
		if(checkEdgeLabel(bi_edges->data[i], bi_edges->data[i+1], u, w, H, G)){
			bi_matrix[j][c_idx] = 1;
		}
		// if we don't have more neighbors in the bi_edges, or next neighbor is the different with current one
		// here it is important to check i+2 is NULL or not, if we change the sequence here then we would access out of the bi_edges
		//if(i+2 >= bi_edges->size || ((struct Vertex*)bi_edges->data[i])->id != ((struct Vertex*)bi_edges->data[i+2])->id){
		if(i+2 < bi_edges->size && ((struct Vertex*)bi_edges->data[i])->number != ((struct Vertex*)bi_edges->data[i+2])->number){
		// increase the row of the bipartite matrix
			++j;
		}
	}
}

/* Print the matrix with specified length of row and column
 **/
void printMatrix(int leng_r, int leng_c, int (*bi_matrix)[leng_c]){
	for(int i=0; i<leng_r; i++){
		for(int j=0; j<leng_c; j++){
			printf(" %i ", bi_matrix[i][j]);
		}
		printf("\n");
	}
}

/* Get the bipartite matching result from here
 **/
void getIsoCharacteristics(int u, struct ShallowGraph* t, struct ShallowGraph* t_prime, Vector* CT, Vector* CTPrime,
		VectorCharacteristics* CC, struct IsoCharacteristics* tmp_iso, struct Vertex* w, int y, struct Graph* H, struct Graph* G){
	Vector neighbors;
	vector_init(&neighbors);
	// here we already remove y if u!=y
	getNeighbor(u, y, &neighbors, H);
	// get the points of edges on the bipartite graph
	Vector bi_edges;
	vector_init(&bi_edges);
	// get the size of bipartite matrix row (number of c)
	int size_theta = CT->size + CTPrime->size;
	// get the vertices of the bipartite graph
	int num_neighbor = get_c_u_primes(&bi_edges, t, t_prime, CT, CTPrime, u, &neighbors, CC);
	// if we have less number of neighbor in Bipartite graph then in neighbor(u) or neighbor(u)/y in H
	// we definitely can not find a match to cover all the neighbors in neighbor(u) or neighbor(u)/y
	if(num_neighbor < neighbors.size){
		vector_free(&neighbors);
		vector_free(&bi_edges);
		return;
	}
	// join c in CT and CTPrime into CTheta
	Vector CTheta;
	vector_init(&CTheta);
	joinCtCtprime(&CTheta, CT, CTPrime);
	// get the bipartite matrix
	int bi_matrix[num_neighbor][size_theta];
	memset(bi_matrix, 0, num_neighbor*size_theta*sizeof bi_matrix[0][0]);
	constructBiMatrix(size_theta, bi_matrix, &bi_edges, &CTheta, u, w, H, G);
	// printMatrix(num_neighbor, size_theta, bi_matrix);
	// construct the returned isoCharacteristics
	tmp_iso->w = w->id;
	// if neighbors match with child of w in t and u same label with w, or we have a single vertex and it has the same label as w
	if((strcmp(H->vertices[u]->label, w->label)==0 && checkMatching(&bi_edges, u, y, w, num_neighbor, size_theta, bi_matrix, H))
			|| (neighbors.size == 0 && strcmp(H->vertices[u]->label, w->label) == 0)){
		tmp_iso->start = y;
		tmp_iso->end = u;
	}
	vector_free(&CTheta);
	vector_free(&bi_edges);
	vector_free(&neighbors);
}

/* Implement the main iterations of characteristics function: We don't combine t and t_prime, but use them together
 * here we keep using the spanning tree generated from the subGraphIsoCheck function (so we can directly check the address
 * of the spanning tree in characteristics of CC with any characteristic generated from this function)
 **/
int isoCharacteristics(struct IsoCharacteristics* tmp_iso, struct Vertex* v, Vector* joi_bi_components, Vector* post_orders, VectorInt* post_order_w, int u, int y,
		struct ShallowGraph* t, struct Vertex* w, VectorCharacteristics* CC, struct Graph* G, struct Graph* H, Vector* spl
){
	struct ShallowGraph* t_prime = NULL;
	Vector CT;
	vector_init(&CT);
	Vector CTPrime;
	vector_init(&CTPrime);
	getChild(&CT, w, post_order_w, G, t);

	// the main iteration of the CHARACTERISTICS function starts from here
	int flag = 0;
	for(int i=0; i<joi_bi_components->size; i++){
		// check whether the vertex w is a rooted vertex in T and not same as v
		if(w->id == ((struct ShallowGraph*)joi_bi_components->data[i])->edges->startPoint->id && w->id != v->id){
			flag = 1;
			// iterate all the spanning trees
			int j = 0;
			for(t_prime=spl->data[i]; t_prime; t_prime=t_prime->next, ++j){
				VectorInt* post_order_w_prime = ((Vector*)(post_orders->data[i]))->data[j];
				// get children of w in graph G
				getChild(&CTPrime, w, post_order_w_prime, G, t_prime);
				getIsoCharacteristics(u, t, t_prime, &CT, &CTPrime, CC, tmp_iso, w, y, H, G);
				if(tmp_iso->start!=-1){
					vector_free(&CT);
					vector_free(&CTPrime);
					return 1;
				}
				// free CTPrime, and re-initialize it to get new data only
				vector_free(&CTPrime);
				vector_init(&CTPrime);
			}
		}
	}
	if(flag == 0){
		getIsoCharacteristics(u, t, t_prime, &CT, &CTPrime, CC, tmp_iso, w, y, H, G);
		if(tmp_iso->start!=-1){
			vector_free(&CT);
			vector_free(&CTPrime);
			return 1;
		}
	}
	vector_free(&CT);
	vector_free(&CTPrime);
	return 0;
}

/* Join two ShallowGraph into one, initially used to join two Bi-connected components,
 * then generating spanning trees from the combined component
 **/
struct ShallowGraph* joinComponent(struct ShallowGraph* g, struct ShallowGraph* h){
	struct VertexList* lhlp = NULL;
	struct VertexList* tail = NULL;
	for(lhlp=g->edges; lhlp; lhlp=lhlp->next){
		tail = lhlp;
	}
	tail->next = h->edges;
	g->m += h->m;
	h->edges = NULL;
	h->lastEdge = NULL;
	return g;
}

/* Get the index of vertex number in the post_order array
 **/
int getIndex(int* post_order, int length, int number){
	int i = 0;
	for(; i<length; i++){
		if(post_order[i] == number){
			break;
		}
	}
	return i;
}

/* Join a list of ShallowGraph which consists of bi-connected components
 **/
void joinBiComponents(Vector* joi_bi_components, struct ShallowGraph* bComponents, struct Graph* G){
	struct ShallowGraph* shlp;
	struct VertexList* e;
	// create an array to hold the shortest distance vertex to r of each 2-connected components
	VectorInt array;
	vint_init(&array);
	// get the post order of the vertices on graph G, with root r=0
	int* post_order = (int*)getPostorder(G, 0);
	int i=0, j;
	// get the shortest distance vertex to r of each 2-connected components into array
	for(shlp=bComponents; shlp; shlp=shlp->next){
		int max = -1;
		for(e=shlp->edges; e; e=e->next){
			int start = getIndex(post_order, G->n, e->startPoint->number);
			if(start > max){
				//max = getIndex(post_order, G->n, e->startPoint->number);
				max = start;
			}
			int end = getIndex(post_order, G->n, e->endPoint->number);
			if(end > max){
				max = end;
			}
		}
		vint_append(&array, post_order[max]);
	}
	//printf("Number of biconnected components: %i\n", array.size);
	// get a copy of array, as an array of root corresponding to each bi component
	VectorInt* array_copy = vint_copy(&array);
	// get distinct vertices in the first NewLength
	int NewLength = vint_distinct(&array);
	// get the combined components in a post order
	for(j=0; j<G->n; j++){
		// iterate all the root vertices
		for(i=0; i<NewLength; i++){
			// if the post order vertex in G is a root vertex
			if(post_order[j] == array.data[i]){
				int count = 0, k = 0;
				// join the 2-connected components with same root v
				for(k=0, shlp=bComponents; k<array.size; k++, shlp=shlp->next){
					if(array_copy->data[k] == array.data[i]){
						if(count == 0){
							vector_append(joi_bi_components, shlp);
						}else{
							joi_bi_components->data[joi_bi_components->size-1] = joinComponent(joi_bi_components->data[joi_bi_components->size-1], shlp);
						}
						++count;
					}
				}
			}
		}
	}
	free(post_order);
	vint_free(&array);
	vint_free(array_copy);
	free(array_copy);
}

/* Print the joined components
 **/
void printJoinComponents(Vector* joi_bi_components){
	struct VertexList *lhlp;
	for(int i=0; i<joi_bi_components->size; i++){
		printf("%i, \n", i);
		for(lhlp=((struct ShallowGraph*)joi_bi_components->data[i])->edges; lhlp; lhlp=lhlp->next){
			printf("%i ", lhlp->startPoint->number);
			printf("%i ", lhlp->endPoint->number);
		}
		printf("\n");
	}
}

/* Get the number of vertex in graph g by their id
 **/
int getNumById(struct Graph* g, int v_id){
	int number = 0;
	for(int i=0; i<g->n; i++){
		if(g->vertices[i]->id == v_id){
			number = g->vertices[i]->number;
			break;
		}
	}
	return number;
}

/* Print elements in structure IsoCharacteristics** array CC
 **/
void printCC(VectorCharacteristics* CC){
	for(int i=0; i<CC->size; i++){
		printCharacteristic(&(CC->data[i]));
		printf("\n");
	}
}

/* Add tmp_iso to current characteristics vector with duplication removed
 **/
void addToChaSet(VectorCharacteristics* CC_current, struct IsoCharacteristics* tmp_iso){
	//if(!vc_contain(CC_current, tmp_iso)){
		vc_append(CC_current, tmp_iso);
	//}
}

/* Calculate the total number of spanning trees from a list of shallow graph cycle
 **/
int numberOfSpanningTrees(Vector* spanning_trees){
	int count = 0;
	for(int i=0; i<spanning_trees->size; i++){
		struct ShallowGraph *idx = spanning_trees->data[i];
		while(idx && idx->next != spanning_trees->data[i]){
			idx = idx->next;
			++count;
		}
		if(idx){
			++count;
		}
		++count;
	}
	return count;
}

void switchPointer(Vector* o, Vector* s){
	s = o;
	o = s->data[0];
}

void assignObjects(struct GraphObjects* object, struct Graph* g, struct ShallowGraph* bComponents, Vector* joi_bi_components, Vector* spanning_tree_list, Vector* sub_graphs, Vector* post_orders){
	object->g = g;
	object->b_components = bComponents;
	object->joi_bi_components = joi_bi_components;
	object->spanning_trees = spanning_tree_list;
	object->sub_graphs = sub_graphs;
	object->post_orders = post_orders;
}

/* Pre processing a single graph
 * 0. get all the 2 connected components in graph
 * 1. get all joined 2 connected components
 * 2. get all locally generated spanning trees
 * 3. get all post orders in each locally generated spanning trees
 */
struct GraphObjects* singleGraphPreProcessing(struct Graph* g, struct ShallowGraph* bComponents, Vector* joi_bi_components, Vector* spanning_tree_list, Vector* sub_graphs, Vector* post_orders, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	struct GraphObjects* object = malloc(sizeof(struct GraphObjects));
	bComponents = listBiconnectedComponents(g, sgp);
	joinBiComponents(joi_bi_components, bComponents, g);

    int count = 0;
	for(int i=0; i<joi_bi_components->size; i++){
		// convert common vertex block to graph
		struct Graph* g = shallowGraphToGraph(joi_bi_components->data[i], gp);
		// get spanning trees from local common vertex graph and get pointer of spanning trees in spanning_tree_list, used globally
		struct ShallowGraph* spanning_tree = listSpanningTrees(g, sgp, gp);
		struct ShallowGraph* shlptree = NULL;
		struct Vertex* v = spanning_tree->edges->startPoint;
		Vector* post_order_in_one_block = malloc(sizeof(Vector));
		vector_init(post_order_in_one_block);
		for(shlptree = spanning_tree; shlptree; shlptree=shlptree->next){
            ++count;
			// get vertex numbers in post order with vertex id
			VectorInt* tmp = getSpanningTreePostOrder(shlptree, gp, v);
			vector_append(post_order_in_one_block, tmp);
		}

		vector_append(post_orders, post_order_in_one_block);
		vector_append(spanning_tree_list, spanning_tree);
		vector_append(sub_graphs, g);
	}
	assignObjects(object, g, bComponents, joi_bi_components, spanning_tree_list, sub_graphs, post_orders);
	return object;
}

/* Post process the list of generated objects used in graph G, used in subtree isomorphic
 * this part frees all objects corresponding to the generated objects in singleGraphPreProcessing()
 */
void singleGraphPostProcessing(struct GraphObjects* objects, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	struct ShallowGraph* bComponents = objects->b_components;
	Vector* joi_bi_components = objects->joi_bi_components;
	Vector* spanning_tree_list = objects->spanning_trees;
	Vector* sub_graphs = objects->sub_graphs;
	Vector* post_orders = objects->post_orders;
	// free locally generated spanning trees
	for(int i=0; i<spanning_tree_list->size; i++){
		dumpShallowGraphCycle(sgp, spanning_tree_list->data[i]);
	}
	vector_free(spanning_tree_list);
	// free generated post_orders
	for(int i=0; i<post_orders->size; i++){
		Vector* post_order_in_one_block = post_orders->data[i];
		for(int j=0; j<post_order_in_one_block->size; j++){
			vint_free(post_order_in_one_block->data[j]);
			free(post_order_in_one_block->data[j]);
		}
		vector_free(post_order_in_one_block);
		free(post_order_in_one_block);
	}
	vector_free(post_orders);
	
	// free all local block graphs
	for(int i=0; i<sub_graphs->size; i++){
		dumpGraph(gp, sub_graphs->data[i]);
	}
	vector_free(sub_graphs);
	
	// free 2 connected components
	dumpShallowGraphCycle(sgp, bComponents);
	// free joi_bi_components
	vector_free(joi_bi_components);
    // free all allocated vectors
	free(joi_bi_components);
	free(spanning_tree_list);
	free(post_orders);
    free(sub_graphs);
    free(objects);
}

/* Subgraph isomorphic check from graph H into graph G
 * All the characteristics from graph H to G can be hold in CC
 **/
int subGraphIsoCheck(struct Graph* H, struct Graph* G, struct GraphObjects* object, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	int X=0;
	// if pattern graph H has more number of vertices, we definitely not be able to find a sub graph isomorphic from graph H to G
	if(H->n > G->n){
		return X;
	}
	int is_single_version = 0;
	struct ShallowGraph* bComponents = NULL;
	Vector* joi_bi_components, *spanning_tree_list, *post_orders;
	// Vector sub_graphs;
    Vector* sub_graphs = NULL;
	// for a single version of subtree isomorphism check
	if(object == NULL){
		printf("Single version\n");
		is_single_version = 1;
		sub_graphs = malloc(sizeof(Vector));
		joi_bi_components = (Vector*)malloc(sizeof(Vector));
		spanning_tree_list = (Vector*)malloc(sizeof(Vector));
		post_orders = (Vector*)malloc(sizeof(Vector));
		vector_init(joi_bi_components); vector_init(spanning_tree_list); vector_init(post_orders); vector_init(sub_graphs);
		object = singleGraphPreProcessing(G, bComponents, joi_bi_components, spanning_tree_list, sub_graphs, post_orders, gp, sgp);
	}else{
		bComponents = object->b_components;
		joi_bi_components = object->joi_bi_components;
		spanning_tree_list = object->spanning_trees;
		post_orders = object->post_orders;
		sub_graphs = object->sub_graphs;
	}
	// hold for characteristics
	VectorCharacteristics CC;
	vc_init(&CC);
	struct ShallowGraph* shlptree = NULL;
	int size = 0;
	// start the iteration of the main function of Algorithm 1, post order of v in V(t)
	for(int i=0; i<joi_bi_components->size; i++){
		size++;
		// get root vertex v of current local block
		struct Vertex* v = ((struct ShallowGraph*)joi_bi_components->data[i])->edges->startPoint;
		Vector* post_order_in_one_block = post_orders->data[i];
		// set of spanning trees rooted at v
		int j = 0;
		for(shlptree=spanning_tree_list->data[i]; shlptree; shlptree=shlptree->next, ++j){
			// get post order of local spanning tree
			VectorInt* post_order_w = post_order_in_one_block->data[j];
			// iterate local spanning tree in post order
			for(int d=0; d<post_order_w->size; d++){
				// get the vertex w in graph G, so we can safely free other temporarily generated graphs
				struct Vertex* w = G->vertices[post_order_w->data[d]];
				// for all the vertices on pattern graph H
				for(int e=0; e<H->n; e++){
					// get the vertex number of e
					int u = H->vertices[e]->number;
					// compute when y == u
					struct IsoCharacteristics* tmp_iso = (struct IsoCharacteristics*)malloc(sizeof(struct IsoCharacteristics));
					initCharacteristic(tmp_iso);
					isoCharacteristics(tmp_iso, v, joi_bi_components, post_orders, post_order_w, u, u, shlptree, w, &CC, G, H, spanning_tree_list);
					// change X=1, as we find a subgraph isomorphic once we get (H(u,u), t, w)
					if(tmp_iso->start != -1){
						// get the post order of current spanning tree with vertex id in graph G
						// tmp_iso->t = shlptree;
						// addToChaSet(&CC, tmp_iso);
						// only print once and return for only subgraph isomorphic check
                        if(is_single_version){
                            singleGraphPostProcessing(object, gp, sgp);
                            printf("Subgraph isomorphic exist!\n");
                        }
						X = 1;
					    free(tmp_iso);
                        vc_free(&CC);
						return X;
					}
					// compute when y are neighbors
					for(struct VertexList* edge=H->vertices[u]->neighborhood; edge!=NULL; edge=edge->next){
						int y = edge->endPoint->number;
						initCharacteristic(tmp_iso);
						isoCharacteristics(tmp_iso, v, joi_bi_components, post_orders, post_order_w, u, y, shlptree, w, &CC, G, H, spanning_tree_list);
						// same as above
						if(tmp_iso->start != -1){
							tmp_iso->t = shlptree;
							addToChaSet(&CC, tmp_iso);
						}
					}
					free(tmp_iso);
				}
			}
		}
	}
	if(is_single_version == 1){
        singleGraphPostProcessing(object, gp, sgp);
	}
	// free characteristics
	vc_free(&CC);
	return X;
}

int subtreeIsomorphic(VectorCharacteristics* CC, struct Graph* H, struct Graph* G, struct GraphObjects* object){
	Vector* joi_bi_components = object->joi_bi_components;
	Vector* spanning_tree_list = object->spanning_trees;
	Vector* post_orders = object->post_orders;
	int X=0;
	// if pattern graph H has more number of vertices, we definitely not be able to find a sub graph isomorphic from graph H to G
	if(H->n > G->n){
		return X;
	}
	struct ShallowGraph* shlptree = NULL;
	int size = 0;
	// start the iteration of the main function of Algorithm 1, post order of v in V(t)
	for(int i=0; i<joi_bi_components->size; i++){
		size++;
		// get root vertex v of current local block
		struct Vertex* v = ((struct ShallowGraph*)joi_bi_components->data[i])->edges->startPoint;
		Vector* post_order_in_one_block = post_orders->data[i];
		// set of spanning trees rooted at v
		int j = 0;
		for(shlptree=spanning_tree_list->data[i]; shlptree; shlptree=shlptree->next, ++j){
			// get post order of local spanning tree
			VectorInt* post_order_w = post_order_in_one_block->data[j];
			// iterate local spanning tree in post order
			for(int d=0; d<post_order_w->size; d++){
				// get the vertex w in graph G, so we can safely free other temporarily generated graphs
				struct Vertex* w = G->vertices[post_order_w->data[d]];
				// for all the vertices on pattern graph H
				for(int e=0; e<H->n; e++){
					// get the vertex number of e
					int u = H->vertices[e]->number;
					// compute when y == u
					struct IsoCharacteristics* tmp_iso = (struct IsoCharacteristics*)malloc(sizeof(struct IsoCharacteristics));
					initCharacteristic(tmp_iso);
					isoCharacteristics(tmp_iso, v, joi_bi_components, post_orders, post_order_w, u, u, shlptree, w, CC, G, H, spanning_tree_list);
					// change X=1, as we find a subgraph isomorphic once we get (H(u,u), t, w)
					if(tmp_iso->start != -1){
						// get the post order of current spanning tree with vertex id in graph G
						tmp_iso->t = shlptree;
						addToChaSet(CC, tmp_iso);
						X = 1;
					}
					// compute when y are neighbors
					for(struct VertexList* edge=H->vertices[u]->neighborhood; edge!=NULL; edge=edge->next){
						int y = edge->endPoint->number;
						initCharacteristic(tmp_iso);
						isoCharacteristics(tmp_iso, v, joi_bi_components, post_orders, post_order_w, u, y, shlptree, w, CC, G, H, spanning_tree_list);
						// same as above
						if(tmp_iso->start != -1){
							tmp_iso->t = shlptree;
							addToChaSet(CC, tmp_iso);
						}
					}
					free(tmp_iso);
				}
			}
		}
	}
	return X;
}

/* Check whether y is between u and a, if it is true, y is above u
 * here we generate post order of Graph rooted at vertex a, check whether y is after u
 **/
int isParent(int u, int a, int y, struct Graph* H){
	// here we get the post order of graph H, rooted at a
	int* post_order = (int*)getPostorder(H, a);
	int u_i=-1, y_i=-1;
	// iterate the H->n number of values in array 'post_order'
	for(int i=0; i<H->n; i++){
		if(post_order[i] == u){
			u_i = i;
		}
		if(post_order[i] == y){
			y_i = i;
		}
	}
	free(post_order);
	// if y is between 'a' and 'u', as here rooted at 'a', so 'a' will certainly be the last one
	if(y_i > u_i){
		return 1;
	}
	return 0;
}

/* Perform fast subgraph mining from graph H into graph G by combing characteristics information from CC_pre
 **/
int fastSubGraphIsoCheck(struct Graph* H, struct Graph* G, struct GraphObjects* object, VectorCharacteristics* CC_pre, VectorCharacteristics* CC_current, int a, int b){
	Vector* joi_bi_components = object->joi_bi_components;
	Vector* spanning_tree_list = object->spanning_trees;
	Vector* post_orders = object->post_orders;
	int X=0;
	// if pattern graph H has more number of vertices, we definitely not be able to find a sub graph isomorphic from graph H to G
	if(H->n > G->n){
		return X;
	}
	struct ShallowGraph* shlptree = NULL;
	// start the iteration of the main function of Algorithm 2, Fast Frequent Subtree Discovery
	for(int i=0; i<joi_bi_components->size; i++){
		// root vertex for current local common vertex block
		struct Vertex* v = ((struct ShallowGraph*)joi_bi_components->data[i])->edges->startPoint;
		// set of spanning trees rooted at v
		int j = 0;
		for(shlptree=spanning_tree_list->data[i]; shlptree; shlptree=shlptree->next, ++j){
			// get post order of local spanning tree
			VectorInt* post_order_w = ((Vector*)post_orders->data[i])->data[j];
			// iterate current spanning tree in a post order fashion
			for(int d=0; d<post_order_w->size; d++){
				// get the vertex of graph G corresponding w in local spanning tree
				struct Vertex* w = G->vertices[post_order_w->data[d]];
				// temporary characteristics used for later all cases
				struct IsoCharacteristics* tmp_iso = (struct IsoCharacteristics*)malloc(sizeof(struct IsoCharacteristics));
				// a == -1, when we have a single vertex
				if(a != -1){
					// Add new characteristics
					if(strcmp(H->vertices[b]->label, w->label) == 0){
						initCharacteristic(tmp_iso);
						constructIsoCharacteristics(tmp_iso, a, b, w, shlptree);
						addToChaSet(CC_current, tmp_iso);
					}
					// check if start_a, end_a in CC_pre
					initCharacteristic(tmp_iso);
					constructIsoCharacteristics(tmp_iso, a, a, w, shlptree);
					if(vc_contain(CC_pre, tmp_iso)){
						initCharacteristic(tmp_iso);
						constructIsoCharacteristics(tmp_iso, b, a, w, shlptree);
						addToChaSet(CC_current, tmp_iso);
					}
				}
				// update CC_current with CC_current UNION CHARACTERISTICS(v,b,b,t,w)
				initCharacteristic(tmp_iso);
				isoCharacteristics(tmp_iso, v, joi_bi_components, post_orders, post_order_w, b, b, shlptree, w, CC_current, G, H, spanning_tree_list);
				if(tmp_iso->start != -1){
					tmp_iso->t = shlptree;
					addToChaSet(CC_current, tmp_iso);
				}
				// if start_b, end_b w in CC_current, set X=TRUE
				initCharacteristic(tmp_iso);
				constructIsoCharacteristics(tmp_iso, b, b, w, shlptree);
				if(vc_contain(CC_current, tmp_iso)){
					// if(X == 0){printf("subgraph isomorphic exist!\n");}
					X=1;
				}
				// Filter existing characteristics w.r.t CC_pre
				if(a != -1){
					for(int idx=0; idx<CC_pre->size; idx++){
						if(CC_pre->data[idx].w == w->id && CC_pre->data[idx].t == shlptree){
							if(isParent(CC_pre->data[idx].end, a, CC_pre->data[idx].start, H)){
								initCharacteristic(tmp_iso);
								constructIsoCharacteristics(tmp_iso, CC_pre->data[idx].start, CC_pre->data[idx].end, w, shlptree);
								addToChaSet(CC_current, tmp_iso);
							}else{
								initCharacteristic(tmp_iso);
								isoCharacteristics(tmp_iso, v, joi_bi_components, post_orders, post_order_w, CC_pre->data[idx].end, CC_pre->data[idx].start, shlptree, w, CC_current, G, H, spanning_tree_list);
								if(tmp_iso->start != -1){
									tmp_iso->t = shlptree;
									addToChaSet(CC_current, tmp_iso);
								}
							}
							// if u u in CC_current, set X=TRUE
							initCharacteristic(tmp_iso);
							constructIsoCharacteristics(tmp_iso, CC_pre->data[idx].end, CC_pre->data[idx].end, w, shlptree);
							if(vc_contain(CC_current, tmp_iso)){
								// if(X == 0){ printf("subgraph isomorphic exist!\n");}
								X=1;
							}
						}
					}
				}
				free(tmp_iso);
			}
		}
	}
	return X;
}

/* Remove vertices from g, which are not in vertices_list */
void removeSingleVertices(struct Graph* g, VectorInt* vl){
	struct Vertex** vertices_copy = malloc(vl->size * sizeof(struct Vertex*));
	for(int i=0; i<g->n; i++){
		for(int j=0; j<vl->size; j++){
			if(vl->data[j] == i){
				vertices_copy[j] = g->vertices[i];
			}
		}
	}
	struct Vertex** tmp = g->vertices;
	g->vertices = vertices_copy;
	free(tmp);
	g->n = vl->size;
}

/* Remove vertex number v in graph g
 **/
void removeSingleVertex(struct Graph* g, struct VertexPool* vp, int v){
	struct Vertex** vertices_copy = malloc((g->n - 1) * sizeof(struct Vertex*));
	for(int i=0, j=0; i<g->n; i++){
		if(i != v){
			vertices_copy[j++] = g->vertices[i];
		}
	}
	dumpVertex(vp, g->vertices[v]);
	struct Vertex** tmp = g->vertices;
	g->vertices = vertices_copy;
	--g->n;
	free(tmp);

	// rearrange vertex number, otherwise some implementation may not work as expected, like canonicalStingOfTree
	for(int i=0; i<g->n; i++){
		g->vertices[i]->number = i;
	}
}

/* Create a new graph g by keeping vertices of graph G in vertices list vl
 **/
struct Graph* createSimpleGraph(struct GraphPool* gp, struct Graph* H, VectorInt* vl){
	struct Graph* g = emptyGraph(H, gp);
	for(int i=0; i<vl->size; i++){
		for(int j=i+1; j<vl->size; j++){
			if(checkNeighbor(H->vertices[vl->data[i]], H->vertices[vl->data[j]])){
				addEdgeBetweenVertices(vl->data[i], vl->data[j], getEdgeLabel(H->vertices[vl->data[i]], H->vertices[vl->data[j]], H), g, gp);
			}
		}
	}
	removeSingleVertices(g, vl);
	return g;
}

/* Create a simple graph from G, by removing vertex number "vertex_remove" from G
 **/
struct Graph* createSGraph(struct GraphPool* gp, struct VertexPool* vp, struct Graph* G, int vertex_remove){
	struct Graph* g = cloneGraph(G, gp);
	struct VertexList* edge = g->vertices[vertex_remove]->neighborhood;
	deleteEdgeBetweenVertices(g, edge, gp);
	removeSingleVertex(g, vp, vertex_remove);
	return g;
	/*
	VectorInt vl;
	vint_init(&vl);
	for(int i=0; i<G->n; i++){
		if(i != vertex_remove){
			vint_append(&vl, i);
		}
	}
	struct Graph* g = createSimpleGraph(gp, G, &vl);
	vint_free(&vl);
	return g;*/
}

/* A singleton implementation of fast subgraph isomorphic checking, from graph H to graph G
 * This is just a way to show the fast version works, in this implementation it is not faster
 * than the naive version, as here we construct the a single vertex graph to graph H, it is not faster than
 * just comparing graph H and graph G
 **/
/*
int fastSubgraphIsoCheck(struct Graph* H, struct Graph* G, Vector* joi_bi_components, Vector* spanning_tree_list, struct GraphPool* gp, struct ShallowGraphPool* sgp){
	int X = 0;
	VectorInt vertices;
	vint_init(&vertices);
	VectorCharacteristics CC_pre;
	vc_init(&CC_pre);
	vint_append(&vertices, H->vertices[0]->number);
	struct Graph* g = createSimpleGraph(gp, H, &vertices);
	if(fastSubGraphIsoCheck(g, G, NULL, &CC_pre, -1, 0 ,gp, sgp)){
		if(H->n == 1){
			X = 1;
		}
	}else{
		return X;
	}
	for(int i=1; i<H->n; i++){
		dumpGraph(gp, g);
		vint_append(&vertices, H->vertices[i]->number);
		g = createSimpleGraph(gp, H, &vertices);
		VectorCharacteristics CC_current;
		vc_init(&CC_current);
		struct Vertex* b = g->vertices[g->n-1];
		struct Vertex* a = b->neighborhood->endPoint;
		if(fastSubGraphIsoCheck(g, G, &CC_pre, &CC_current, a->number, b->number, gp, sgp)){
			if(g->n == H->n){
				X = 1;
				printf("Subgraph Isomorphic Exist!\n");
			}
		}
		vc_free(&CC_pre);
		CC_pre = CC_current;
		printCC(&CC_current);
	}
	dumpGraph(gp, g);
	vc_free(&CC_pre);
	vint_free(&vertices);
	return X;
}
*/
