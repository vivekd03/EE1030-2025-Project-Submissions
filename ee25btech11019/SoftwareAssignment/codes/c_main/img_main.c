#include<stdio.h>
#include<stdlib.h>
#include"img_svd.h"

int main(){
    int m, n;//height,width
    const char* filename = "inputbhai.jpg";
    double **A = read_image_grayscale(filename, &n, &m);

    if(A == NULL){//jo jpg nathi to png try kariye
        filename = "inputbhai.png";
        printf("Couldn't find .jpg, trying .png\n");
        A = read_image_grayscale(filename, &n, &m);
    }

    if(A == NULL){
        fprintf(stderr,"Failed to read input image.\n");
        return 1;
    }

    int k;//input k lyo
    printf("Enter the number of singular values (k): ");
    scanf("%d",&k);
    int k_max = (n < m) ? n : m;//k max rank thi toh nano j hovo joiye max rnak = min(dimension)
    if(k>k_max){
        printf("Changing k (%d) to max possible (%d).\n", k, k_max);
        k =k_max;
    }
    if(k <=0){
        printf("Error: k must be positive.\n");
        return 1;
    }

    printf("Allocating memory...\n");
    double **U =creat_mat(m, k);// (m k)
    double **V =creat_mat(k, n);// (k n) (V_kt)
    double *sing =malloc(k*sizeof(double)); // k x 1
    double **A_k =creat_mat(m, n);//(m n)(reconstructed matrix)

    printf("Calculating top %d singular values \n",k);
    topkSingular(A, m,n, k, U, V, sing);
    
    for(int i =0;i<m; i++){//A_k[i][j] = sum of (l =0; l<k) (U[i][l]*sing[l]*V[l][j]) A_k ma U and V nakhva
        for(int j =0;j<n; j++){
            double sum =0.0;
            for(int l =0;l<k; l++){
                sum += U[i][l]*sing[l]*V[l][j];
            }
            A_k[i][j] =sum;
        }
    }
    
    printf("\nimage done\n");
    write_image_grayscale("output_test.jpg", A_k, n,m, 90);//90 quality A_k ma changes karya che A ma nai
    
    printf("Frobenius norm: %lf\n", frobn(A, A_k, m, n));

    printf("freeing memory.\n");
    free_mat(A, m),free_mat(A_k, m),free_mat(U, m),free_mat(V, k),free(sing);
    printf("Done.\n");
    return 0;
}