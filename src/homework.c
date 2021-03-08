
#include"fem.h"



# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->mesh  = femMeshRead(filename);
    theProblem->edges = femEdgesCreate(theProblem->mesh);
    if (theProblem->mesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theProblem->mesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theProblem->mesh->nNode);
    return theProblem;
}

# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    free(theProblem);
}


# endif
# ifndef NOMESHLOCAL


void femMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y)
{
    if (theMesh->nLocalNode == 3){
        x[0] = theMesh->X[theMesh->elem[i*3]] ;
        x[1] = theMesh->X[theMesh->elem[i*3+1]] ;
        x[2] = theMesh->X[theMesh->elem[i*3+2]] ;

        y[0] = theMesh->Y[theMesh->elem[i*3]] ;
        y[1] = theMesh->Y[theMesh->elem[i*3+1]] ;
        y[2] = theMesh->Y[theMesh->elem[i*3+2]] ;

        map[0] = theMesh->elem[i*3] ;
        map[1] = theMesh->elem[i*3+1] ;
        map[2] = theMesh->elem[i*3+2] ;
    }

    if (theMesh->nLocalNode == 4){
        x[0] = theMesh->X[theMesh->elem[i*4]] ;
        x[1] = theMesh->X[theMesh->elem[i*4+1]] ;
        x[2] = theMesh->X[theMesh->elem[i*4+2]] ;
        x[3] = theMesh->X[theMesh->elem[i*4+3]] ;

        y[0] = theMesh->Y[theMesh->elem[i*4]] ;
        y[1] = theMesh->Y[theMesh->elem[i*4+1]] ;
        y[2] = theMesh->Y[theMesh->elem[i*4+2]] ;
        y[3] = theMesh->Y[theMesh->elem[i*4+3]] ;

        map[0] = theMesh->elem[i*4] ;
        map[1] = theMesh->elem[i*4+1] ;
        map[2] = theMesh->elem[i*4+2] ;
        map[3] = theMesh->elem[i*4+3] ;
    }

}

# endif
# ifndef NOPOISSONSOLVE


void femPoissonSolve(femPoissonProblem *theProblem)
{

    femMesh *mesh = theProblem->mesh;
    femFullSystem *system = theProblem->system ;
    femEdges *edges = theProblem->edges;
    femDiscrete *spacex = theProblem->space;
    femIntegration *rules = theProblem->rule;

    double **A = system->A;
    double *B = system->B;

    double *X = mesh->X;
    double *Y = mesh->Y ;

    int node = mesh->nLocalNode ;

    double x[node];
    double y[node];
    int map[node];

    double dxsi = 0.0 ;
    double dxeta = 0.0 ;
    double dyxsi = 0.0 ;
    double dyeta = 0.0 ;

    double J = 0.0;

    double phi[node] ;
/*
    double xsi = 0 ;
    double eta = 0 ;
*/
    double dphidxsi[node] ;
    double dphideta[node] ;

    double dphidx[node] ;
    double dphidy[node] ;

    for (int i = 0 ; i < mesh->nElem ; i++) {
        femMeshLocal(mesh, i, map, x, y);
        for (int ii = 0 ; ii < rules->n ; ii++) {
            femDiscretePhi2(spacex,rules->xsi[ii],rules->eta[ii],phi);
            femDiscreteDphi2(spacex, rules->xsi[ii], rules->eta[ii], dphidxsi, dphideta);

            //Reset variables
            dxsi = 0.0;
            dxeta = 0.0;
            dyxsi = 0.0;
            dyeta = 0.0;
            //xsi = 0 ;
            //eta = 0 ;


            for (int j = 0; j < node; j++) {

                dxsi += x[j] * dphidxsi[j];
                dyxsi += y[j] * dphidxsi[j];
                dxeta += x[j] * dphideta[j];
                dyeta += y[j] * dphideta[j];
            }

            J = fabs(dxsi * dyeta - dxeta * dyxsi);
            //printf("JAcobien : %f avec dxsi : %f,dyxsi : %f, dxeta : %f, dyeta : %f \n",J, dxsi, dyxsi, dxeta,dyeta) ;

            for (int k = 0; k < spacex->n; k++) {
                dphidx[k] = (1.0 / J) * (dyeta * dphidxsi[k] - dyxsi * dphideta[k]);
                dphidy[k] = (1.0 / J) * (dxsi * dphideta[k] - dxeta * dphidxsi[k]);
            }

            for (int li = 0; li < spacex->n; li++) {
                for (int co = 0; co < spacex->n; co++) {
                    /*
                    double q = rules->weight[li]*J*(dphidx[li]*dphidx[co]+dphidy[li]*dphidy[co]) ;
                    printf("%f\n",q) ;
                     */
                    A[map[li]][map[co]] += rules->weight[li] * J * (dphidx[li] * dphidx[co] + dphidy[li] * dphidy[co]);

                }
                B[map[li]] += J * phi[li] * rules->weight[li] ;
            }
        }
    }

    for (int i= 0; i < edges->nEdge; i++) {
        if (edges->edges[i].elem[1] == -1) {
            int n1 = edges->edges[i].node[0];
            int n2 = edges->edges[i].node[1]; //node[1] ??
            femFullSystemConstrain(system,n1,0);
            femFullSystemConstrain(system,n2,0);
        }
    }

    //femFullSystemPrint(system) ;

    femFullSystemEliminate(system) ;
}

# endif
