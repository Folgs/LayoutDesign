#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 9
#define P 3 /*N number of points, P-1 dimension to draw*/


double modul( double [] );

/* To compute de top eigenvectors of D^-1 A.
INPUT:  in "entrada.in", the weighted adjacendy matrix.
OUTPUT: in "sortida.out", P-1 coumns, the eigenvectors.
NEEDED: You need Couenne and AMPL.*/

int main(void) {

   int l, i, j, iter=0;
   double u1[P-1][N], u2[P-1][N], u[P][N], mod, num, ini[N];
   double tol = 1-1.e-17;
   double prod[2];
   double D[N][N], w[N][N];
   char ent[]="entrada.in", sor[]="sortida.out", cou[]="ortonormal.mod", csor[]="proau.out";

   FILE *entrada, *sortida, *orton, *cout;

   entrada = fopen(ent, "r");
   sortida = fopen(sor, "w");

   for (j = 0; j < 2; j++) {
      for (i = 0; i < N; i++)
         u[j][i] = 1.;
      mod = modul(u[j]);
      for (i = 0; i < N; i++)
         u[j][i] = u[j][i] / mod;
   }

   /*Read weighted adjacency matrix*/
	for (i = 0; i < N; i++) {
		D[i][i] = 0;
      for (j = 0; j < N; j++) {
			fscanf(entrada, "%le ", &w[i][j]);
         D[i][i] = D[i][i] + w[i][j];
		}
	}

   /*Check of matrix (symmetric?)*/
   for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
         if (w[i][j] != w[j][i]) {
            printf("Weight matrix may not symmetric\n");
            exit(1);
         }
      }
   }

   for (j = 0; j < N; j++) {
      fscanf(entrada, "%le ", &ini[i]);
      u2[0][j] = ini[j];
      u2[1][j] = ini[j];
   }

   /*D-orthonormalize*/
   orton = fopen(cou, "w");
      /*Variables creation*/
   for (i = 0; i < N; i++)
      fprintf(orton, "var x%d:=%le;\nvar y%d:=%le;\n", i+1, u2[0][i], i+1, u2[1][i]);
      /*Minimitzatiion*/
   fprintf(orton, "minimize obj: 0");
   for (i = 0; i < N; i++)
      fprintf(orton, "+ (x%d - %le)^2 + (y%d - %le)^2", i+1, u2[0][i], i+1, u2[1][i]);
      /*conditions*/
   fprintf(orton, ";\nsubject to c1: 0");
   for (i = 0; i < N; i++)
      fprintf(orton, "+%le*x%d", D[i][i], i+1);
   fprintf(orton, "=0;\nsubject to c2: 0");
   for (i = 0; i < N; i++)
      fprintf(orton, "+%le*y%d", D[i][i], i+1);
   fprintf(orton, "=0;\nsubject to c3: 0");
   for (i = 0; i < N; i++)
      fprintf(orton, "+%le*x%d*y%d", D[i][i], i+1, i+1);
   fprintf(orton, "=0;\n");
      /*Pass to Couenne*/
   fclose(orton);
   system("ampl <couen.run");
      /*Save the data*/
   cout=fopen(csor,"r");
   for (i = 0; i < N; i++)
      fscanf(cout, "%le %le", &u2[0][i], &u2[1][i]);

   /*Normalization*/
   for (i = 0; i < (P-1); i++) {
      mod = modul(u2[i]);
      for (j = 0; j < N; j++)
         u2[i][j] = u2[i][j] / mod;
   }

   prod[0] = 0.; prod[1] = 0.; iter = 0;

   while (iter < 1000 && (prod[0] < tol || prod[1] < tol)) {
      printf("iter=%d prod-1 =(%le,%le)\n", iter, prod[0]-1, prod[1]-1);
      
      iter = iter + 1;

      for (i = 0; i < N; i++) {
         for (j = 0; j < (P-1); j++)
	    u1[j][i] = u2[j][i];
      }

      /*D-orthonormalize*/
      orton = fopen(cou, "w");
         /*Variable creation*/
      for (i = 0; i < N; i++)
         fprintf(orton, "var x%d:=%le;\nvar y%d:=%le;\n", i+1, u1[0][i], i+1, u1[1][i]);
         /*Minimitzation*/
      fprintf(orton, "minimize obj: 0");
      for (i = 0; i < N; i++)
         fprintf(orton, "+ (x%d - %le)^2 + (y%d - %le)^2", i+1, u1[0][i], i+1, u1[1][i]);
         /*conditions*/
      fprintf(orton, ";\nsubject to c1: 0");
      for (i = 0; i < N; i++)
         fprintf(orton, "+%le*x%d", D[i][i], i+1);
      fprintf(orton, "=0;\nsubject to c2: 0");
      for (i = 0; i < N; i++)
         fprintf(orton, "+%le*y%d", D[i][i], i+1);
      fprintf(orton, "=0;\nsubject to c3: 0");
      for (i = 0; i < N; i++)
         fprintf(orton, "+%le*x%d*y%d", D[i][i], i+1, i+1);
      fprintf(orton, "=0;");
        /* Pass to Couenne*/
      fclose(orton);
      system("ampl <couen.run");
         /*Save the data*/
      cout = fopen(csor, "r");
      for (i = 0; i < N; i++)
         fscanf(cout, "%le %le", &u1[0][i], &u1[1][i]);
      fclose(cout);

      /*Multiply with 0.5(I+D-1A)*/
      for (l = 0; l < (P-1); l++) {
         for (i = 0; i < N; i++) {
            num = 0;
            for (j = 0; j < N; j++) {
               if (w[i][j] != 0)
                  num = num + w[i][j]*u1[l][j];
            }
            u2[l][i] = 0.5*(u1[l][i] + num/D[i][i]);
         }
      }

      /*Normalization*/
      for (j = 0; j < (P-1); j++) {
         mod = modul(u2[j]);
         for (i = 0; i < N; i++)
            u2[j][i] = u2[j][i] / mod;
      }

      /*Compute direction variation*/
      prod[0] = 0; prod[1] = 0;
      for (j = 0; j < (P-1); j++) {
         for (i = 0; i < N; i++)
            prod[j] = prod[j] + u1[j][i]*u2[j][i];
      }

      /*Save the resultats*/
      for (j = 0; j < (P-1); j++) {
         for (i = 0; i < N; i++)
            u[j+1][i] = u1[j][i];
      }

	}


   printf("Final amb %d iteracions\n",iter);

   /*Write the final results*/
   for (i = 0; i < N; i++) {
      for (j = 1; j < P; j++)
         fprintf(sortida, "%le ", u[j][i]);
      fprintf(sortida, "\n");
   }

   printf("prodFinal=(%le,%le)\n", prod[0], prod[1]);

   return 0;
}

/*Compute the modulus of vector v of lenght N
INPUT: A vector of lenght N
OUTPUT: Its modulus*/
double modul( double v[N] ) {
   double mod = 0;
   int i;

   for (i = 0; i < N; i++)
      mod = mod + v[i]*v[i];

   return sqrt(mod);
}
