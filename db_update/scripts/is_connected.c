#include <igraph/igraph.h>

int main(int argc, char **argv) {
    igraph_bool_t res;
    igraph_t g;
    FILE *input;

    /* Without names and weights */
    input=fopen(argv[1], "r");
    if (!input) {
        return 1;
    }

    igraph_read_graph_ncol(&g, input, NULL, 0, IGRAPH_ADD_WEIGHTS_YES, 0);
    fclose(input);

    long int v = igraph_vcount(&g);
    long int e = igraph_ecount(&g);

    igraph_is_connected(&g, &res, 2);
    printf("%ld\t%ld\t%s\n", v, e, res ? "TRUE" : "FALSE");
    igraph_destroy(&g);
    return 0;
}
