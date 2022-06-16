#include <igraph/igraph.h>
#include <math.h>

int compare (const void * a, const void * b)
{
        if (*(double*)a > *(double*)b) return 1;
        else if (*(double*)a < *(double*)b) return -1;
        else return 0;
}


igraph_vector_t find_index(igraph_vector_t *v1, igraph_vector_t *v2, double *s){
        static long int i;
        static double k;
        for (i=0; i<igraph_vector_size(v1); i++) {
                k = (double) VECTOR(*v1)[i];
                if (fabs(k) <= fabs(*s)) {
                        //printf("%f %f\n", *s, k);
                        igraph_vector_push_back(v2, i);
                }
        }
        return *v2;
}


void print_vector(igraph_vector_t *v, FILE *f) {
        long int i;
        for (i=0; i<igraph_vector_size(v); i++) {
                fprintf(f, " %lf", (double) VECTOR(*v)[i]);
        }
        fprintf(f, "\n");
}


igraph_t get_sub(igraph_t *g, igraph_vector_t *wu, igraph_vector_t *wsort, long int *i){
        static igraph_t g1;
        static igraph_vector_t indices;
        igraph_vector_init(&indices, 0);
        static igraph_es_t es;
        igraph_copy(&g1, g);
        static double w;
        w = (double) VECTOR(*wsort)[*i];
        find_index(wu, &indices, &w);
        es = igraph_ess_vector(&indices);
        igraph_delete_edges(&g1, es);
        igraph_vector_destroy(&indices);
        return g1;
}


long int binarySearch(igraph_t *g, igraph_vector_t *v, igraph_vector_t *vsort,  long int l,  long int r){
        static igraph_bool_t res;
        static igraph_t G;
        //igraph_copy(&g2,G1);
        static long int mid;


        if (l > r)
        {
                return -1;
        }

        mid = (l+r)/2;
        G = get_sub(g, v, vsort, &mid);
        igraph_is_connected(&G, &res, 2);
        igraph_integer_t c;
        igraph_clusters(&G, NULL, NULL, &c, IGRAPH_STRONG);

        //printf("%s\tmin: %ld\tmid: %ld\tmax: %ld\tmid weight: %f\n", "STEP0", l, mid, r,  igraph_vector_e(vsort, mid));
        if ((mid - l == 0) && res) {
                //m1 = mid;
                //printf("%s\tmin: %ld\tmid: %ld\tmax: %ld\tmid weight: %f\n", "STEP0", l, mid, r,  igraph_vector_e(vsort, mid));
                //G = get_sub(g, v, vsort, &m);
                return mid;
        }else if ((r - l == 0) && res) {
                //m1 = mid - 1;
                //printf("%s\tmin: %ld\tmid: %ld\tmax: %ld\tmid weight: %f\n", "STEP1", l, mid, r,  igraph_vector_e(vsort, mid));
                //G = get_sub(g, v, vsort, &m);
                return mid;
        }else if (mid != l && !res) {
                //printf("%s\tmid: %ld\tmin: %ld\tmax: %ld\n", "STEP3", mid, l, r);
                //printf("%s\tmin: %ld\tmid: %ld\tmax: %ld\tmid weight: %f\tcomp: %d\n", "STEP3", l, mid, r,  igraph_vector_e(vsort, mid), c);
                igraph_vector_t v1, vsort1;
                igraph_vector_copy(&vsort1, vsort);
                igraph_vector_copy(&v1, v);
                igraph_destroy(&G);
                binarySearch(g, &v1, &vsort1, l, mid);
        }else if (mid != l && res) {
                //printf("%s\tmid: %ld\tmin: %ld\tmax: %ld\n", "STEP4", mid, l, r);
                //printf("%s\tmin: %ld\tmid: %ld\tmax: %ld\tmid weight: %f\tcomp: %d\n", "STEP4", l, mid, r,  igraph_vector_e(vsort, mid), c);
                igraph_vector_t v1, vsort1;
                igraph_vector_copy(&vsort1, vsort);
                igraph_vector_copy(&v1, v);
                igraph_destroy(&G);
                binarySearch(g, v, vsort, mid, r);
        }else{
                return -1;
                igraph_destroy(&G);
        }
        igraph_destroy(&G);
        return mid;
}

int main(int argc, char **argv) {
        igraph_i_set_attribute_table(&igraph_cattribute_table);
        long int nedges, nertex, m, maxw, minw, max;
        int n, k = 0;
        const char *weight = "weight";
        const char *name = "name";
        igraph_t g, g1;
        igraph_vector_t weightatt, weightattsort, weightattsort_u, edges, evc;
        igraph_strvector_t names;
        igraph_real_t *weightatt_c, d, value;
        igraph_integer_t c;
        double w;
        igraph_bool_t res;
        igraph_es_t es;
        long int v,e;
        igraph_arpack_options_t options;
        igraph_arpack_options_init(&options);
        char *nam;
        const char * format = "\n# Results:\n\nCut_position\t%lu\nCut_weight\t%f\nMin_weight\t%f\nMax_weight\t%f\nNum_vrtx\t%ld\nNum_edges\t%ld\nDensity\t%f\nRepresentative\t%s\nEVcent_rep\t%f\nNum_components\t%d\nConnected\t%s\n";

        static igraph_vector_t indices;
        printf("%s\n", "Reading graph...");
        FILE *input;
        /* Without names and weights */
        input = fopen(argv[1], "r");
        if (!input) {
                return 1;
        }
        igraph_read_graph_ncol(&g, input, NULL, 1, 1, 0);
        fclose(input);

        printf("%s\n", "Copying graph...");
        igraph_copy(&g1, &g);

        nedges=igraph_ecount(&g);
        nertex=igraph_vcount(&g);

        printf("%s\n", "Initializing egde vectors...");
        /*Declare vectors and initialise*/
        igraph_vector_init(&weightatt,0);
        igraph_vector_init(&weightattsort,0);

        printf("%s\n", "Getting egde attributes...");
        /* Get weights */
        igraph_cattribute_EANV(&g,weight,igraph_ess_all(IGRAPH_EDGEORDER_ID), &weightatt);

        printf("%s\n", "Converting igraph vector to C array...");
        weightatt_c=(igraph_real_t*) malloc(nedges* sizeof(igraph_real_t));
        igraph_vector_copy_to(&weightatt, weightatt_c);

        printf("%s\n", "Sorting egde vectors...");
        qsort(weightatt_c, nedges, sizeof(double), compare);

        printf("%s\n", "Converting sorted C array to igraph vector...");
        igraph_vector_view(&weightattsort, weightatt_c, nedges);

        printf("%s\n", "Removing duplicates from vector...");
        for( n = 0; n < nedges; n++ )
        {
                if(k==0)
                        weightatt_c[k++]=weightatt_c[n];
                else
                {
                        if(weightatt_c[n]==weightatt_c[k-1])
                                continue;
                        else
                                weightatt_c[k++]=weightatt_c[n];
                }
        }


        printf("%s\n", "Converting unique sorted C array to igraph vector...");
        igraph_vector_view(&weightattsort_u, weightatt_c, nedges);

        printf("%s\n", "Finding weight where to filter...");
        m = binarySearch(&g, &weightatt, &weightattsort_u, 0, igraph_vector_size(&weightattsort_u));

        printf("Selected weight position: %ld\n", m);

        if (m == -1) {
                printf("%s\n", "ERROR: Couldn't find a weight to filter");
                exit(1);
        }

        printf("%s\n", "Removing edges...");
        igraph_vector_init(&indices, 0);
        igraph_vector_init(&edges,0);
        w = (double) VECTOR(weightattsort_u)[m];
        find_index(&weightatt, &indices, &w);
        igraph_es_vector(&es, &indices);
        igraph_delete_edges(&g, es);
        igraph_is_connected(&g, &res, 2);

        if (!res) {
          printf("%s\n", "Resulting graph not connected. Restoring original graph...");
          igraph_destroy(&g);
          igraph_copy(&g, &g1);
          igraph_is_connected(&g, &res, 2);
        }

        printf("%s\n", "Calculating eigenvector centrality...");
        igraph_strvector_init(&names, 0);
        igraph_vector_destroy(&weightatt);
        igraph_vector_init(&weightatt,0);
        printf("%s\n", "HEER...");
        igraph_cattribute_EANV(&g,weight,igraph_ess_all(IGRAPH_EDGEORDER_ID), &weightatt);
        igraph_cattribute_VASV(&g, name, igraph_vss_all(), &names);
        igraph_vector_init(&evc, 0);
        igraph_eigenvector_centrality(&g, &evc, &value, /*directed=*/ 0,
                                      /*scale=*/ 1, &weightatt,  &options);
        max = igraph_vector_which_max(&evc);
        igraph_strvector_get(&names, max, &nam);

        printf("%s\n", "Writing trimmed graph...");
        FILE *output;
        output=fopen("trimmed_graph.ncol", "w");
        igraph_write_graph_ncol(&g, output, "name", "weight");
        //igraph_write_graph_graphml(&g1, output, 0);
        fclose(output);

        printf("%s\n", "Calculating some statistics...");
        igraph_density(&g, &d, 0);
        maxw = igraph_vector_which_max(&weightatt);
        minw = igraph_vector_which_min(&weightatt);
        v = igraph_vcount(&g);
        e = igraph_ecount(&g);
        igraph_clusters(&g, NULL, NULL, &c, IGRAPH_STRONG);
        printf(format, m, w, igraph_vector_e(&weightatt, minw), igraph_vector_e(&weightatt, maxw),
               v, e, d, nam, igraph_vector_e(&evc, max), c, res ? "TRUE" : "FALSE");

        igraph_strvector_destroy(&names);
        igraph_vector_destroy(&evc);
        igraph_vector_destroy(&weightatt);
        igraph_vector_destroy(&weightattsort);
        igraph_destroy(&g);
        igraph_destroy(&g1);
        return 0;
}

