# master-thesis

##### Description
This repository contains the source code of thesis project implementations.

It is part of another framework so you cannot compile it only with these codes.

The compiled program [subgraph](https://github.com/liyakun/master-thesis/blob/master/subgraph) is provided.

##### Help
`subgraph -h`

##### Example

###### Subtree Isomorphism Check
`./subgraph -u isomorphic ./data/isomorphic_example.txt`

###### Frequent Connected Subtree Mining

`./subgraph -u mining -m levelwise -i iterative -s 1 -v 5 -c prefix ./data/mining_example.txt`
