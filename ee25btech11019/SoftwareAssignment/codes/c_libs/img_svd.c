#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

double **creat_mat(int m, int n){
    double** M = (double **)malloc(m*sizeof(double *));// aaiya size of double* aave 2D array ma 1D array rakhvano che
    if(M == NULL){// jo malloc ma jagya no hoy to NULL aape
        exit(1);// no matlab atyarej abnormal rite terminate karo code ne and exit(0) normal jem haltu hoy em halva do
    }
    for(int i = 0; i < m; i++){
        M[i] = (double *)malloc(n*sizeof(double));
        if(M[i] == NULL){
            exit(1);
        }
    }
    return M;
}

void free_mat(double ** M, int m){
    for(int i =0 ;i< m; i++){
        free(M[i]);
    }
    free(M);
}

// C = AB where A (m p)  B (p n)  C (m n)

void mult_mat(double **A, double **B, double** C, int m, int p, int n){
    for(int i=0; i<m; i++){
        for(int j =0; j<n; j++){
            double sum = 0.00;
            for(int k=0;k<p; k++){
                sum+= A[i][k]*B[k][j];
            }
            C[i][j] = sum;
        }
    }
}

void trans_mat(double **A, double **At, int m, int n){
    for(int i =0; i<m; i++){
        for(int j =0; j<n; j++){
            At[j][i] = A[i][j];
        }
    }
}

double vec_dot(double *a, double*b, int n){
    double sum =0.0;
    for(int i =0; i<n; i++){
        sum += a[i]*b[i];
    }
    return sum;
}

void vec_devide(double *a, int n, double k){
    for(int i =0; i<n; i++){
        a[i] /= k;
    }
}

void vecmatMult(double **A, int m, int n , double*x, double *y){// A(m n) x [vector] (n 1) [y = Ax];
    for(int i = 0; i <m; i++){
        double sum = 0.0; 
        for(int j =0; j<n; j++){
            sum += A[i][j]*x[j];
        }
        y[i] = sum;
    }
}

void gramSchmidt(double *v,double **old_V, int n, int count){//count no matlab ketla orthonormal vector che; old_V[0]...old_V[count-1];
    if(count == 0){
        return;
    }
    
    for(int i = 0; i<count; i++){// oethogonalize against all previous vectors
        double p = vec_dot(old_V[i], v,n);//projection lenght
        for(int j =0; j<n; j++){
            v[j] = v[j] - p*old_V[i][j];//jo nahi hoy to p =0;
        }
    }


    double norm = sqrt(vec_dot(v,v,n));// normalize
    if(norm > 1e-12){
        vec_devide(v, n, norm);
    }else{
        // jo vector zero bani jay orthogonalization pachi to ek nani value set kari daiye
        for(int j = 0; j<n; j++){
            v[j] = 1e-10;
        }
        norm = sqrt(vec_dot(v, v,n));
        vec_devide(v,n, norm);
    }
}

void power_iteration(double**A, int m, int n, double *v, double *sing, double** old_V, int count, int maxI, double tolerence){
    double *Av = malloc(m*sizeof(double));
    double *z = malloc(n*sizeof(double));// z = (AtAv)
    for(int i =0; i<n; i++){
        v[i] = 1.0;
    }
    double norm = sqrt(vec_dot(v, v,n));
    vec_devide(v, n, norm);// normalize 
    double old_sing = 0.0;

    for(int iter =0; iter <maxI ; iter++){//i ni jagya a iter karyu!!!!!!
        vecmatMult(A, m,n, v, Av);// y = Av;
        // double At[n][m];
        // trans_mat(A,At,m,n);
        // vecmatmult(At, m, n, Av, z);// z = Aty = (AtA)v;
        for(int j = 0;j<n;j++){
            double sum = 0.0;
            for(int i=0;i<m;i++){
                sum += A[i][j]*(Av[i]);// i vapryu row mate je pela k hatu!!!!!!
            }
            z[j] = sum;
        }// uper comment karyu e same kaam without extra memory

        // z badha juna thi orthogonal hoy
        gramSchmidt(z, old_V, n, count);

        double znorm = sqrt(vec_dot(z, z, n));
        if(znorm ==0.0) break;
        vec_devide(z, n, znorm);// normalize z;

        double maxdiff = 0.0;
        for(int i = 0;i <n; i++){
            double d = fabs(z[i] - v[i]);
            if(d > maxdiff){
                maxdiff = d;
            }
        }

        // v ne update karvo padse aagad na ma y = Av karva mate v k +1 valu yaad kar
        for(int i =0; i<n; i++){
            v[i] =  z[i];
        }

        if(maxdiff<tolerence){//convergence check
            break;
        }
    }
    free(Av),free(z);
}

void topkSingular(double **A, int m, int n, int k, double **U, double **V, double *sing){//computing top k singular values
    if(k <= 0) return;
    if(k >n) k = n;

    double **A_copy = creat_mat(m,n);// A ni copy banavi lidhi
    for(int i =0; i< m; i++){
        for(int j =0;j<n; j++){
            A_copy[i][j] = A[i][j];
        }
    }
       
    // U, V, sing are now passed in from main
    double *v = malloc(n*sizeof(double));
    double *Av = malloc(m*sizeof(double));

    for(int i = 0; i<k; i++){// k singular values joiye che
        power_iteration(A_copy, m, n, v, &sing[i], V, i, 2000, 1e-9);// copy vaprvi

        vecmatMult(A_copy, m,n, v, Av);// sigma = | A*v |;
        double s2 = vec_dot(Av, Av, m);
        double sigma;
        if(s2> 0.0){
            sigma = sqrt(s2);
        }else{
            sigma = 0.0;
        }
        sing[i] = sigma;

        if(sigma > 1e-12){// u (A*v) /sigma;
            vec_devide(Av, m, sigma);
        }else{
            for(int j =0;j<m; j++){ 
                Av[j] = 0.0;
            }
        }
        
        //u & v ne U & V ma nakhva
        for(int j =0;j<n; j++){
            V[i][j] = v[j];// V(k n) - This is V_k^T
        }
        for(int j =0;j<m; j++){
            U[j][i] = Av[j];// U(m k)
        }

        if(sigma >1e-12){
            for(int j =0; j<m; j++){//Acopy ne deflate karo
                for(int l =0;l<n; l++){// k thi confusion no thay
                    A_copy[j][l] = A_copy[j][l] - sing[i]*U[j][i]*V[i][l];//deflate
                }
            }
        }
    }

    // All printf and free calls removed. Main will handle this.
    free_mat(A_copy, m),free(v),free(Av);
}

double** read_image_grayscale(const char *filename, int *width, int *height){
    int n_channels;
    unsigned char *img_data = stbi_load(filename, width, height, &n_channels, 1);//n channel hoy pan aapde forcefully 1 karie(greyscale)
    
    if(img_data == NULL){
        fprintf(stderr,"Error: Failed to load image %s\n",filename);
        return NULL;//double ptr che 1 na kari saki
    }

    printf("Successfully read %s (%dx%d).\n", filename, *width, *height);

    double** matrix =creat_mat(*height,*width);

    for(int y =0;y< *height; y++){//unsigned char buffer thi double matrix;
        for(int x =0;x< *width; x++){
            matrix[y][x] =(double)img_data[y*(*width) + x];
        }
    }
    stbi_image_free(img_data);//stbi_image na chale
    return matrix;
}

void write_image_grayscale(const char *filename, double **matrix, int width,int height, int quality){
    // double min_val =matrix[0][0];
    // double max_val =matrix[0][0];
    // for(int y =0;y<height; y++){
    //     for(int x =0;x<width; x++){
    //         if(matrix[y][x]<min_val){
    //             min_val =matrix[y][x];
    //         }
    //         if(matrix[y][x]>max_val){
    //             max_val =matrix[y][x];
    //         }
    //     }
    // }

    // double range = max_val - min_val;
    // if(range<1e-9){//flat matrix handle kare
    //     range =1.0; 
    // }

    unsigned char *img_data =(unsigned char*)malloc(width*height*sizeof(unsigned char));//buffer allocate karie chi aakhi image mate
    if(!img_data){// == NULL ptr na hale
        fprintf(stderr,"Error: Malloc failed for image buffer.\n");
        return;
    }

    for(int y =0;y<height; y++){//buffer ne fill
        for(int x =0;x<width; x++){
            double norm_val = matrix[y][x];
            if(norm_val <0.0){
                norm_val =0.0;
            }
            if(norm_val >255.0){
                norm_val =255.0;
            }
            img_data[y*width + x] =(unsigned char)(norm_val);
        }
    }

    if(stbi_write_jpg(filename, width, height,1, img_data, quality) == 0){// l channels write image
        fprintf(stderr,"Error: Failed to write image.\n");
    }

    free(img_data);//image buffer free
}

double frobn(double **A,double **B, int m,int n){//khali 2 norm j che pan matrix mate
    double sum =0.0;
    for(int i =0;i<m; i++){
        for(int j =0;j<n; j++){
            double f = (A[i][j]-B[i][j]);
            sum += f*f;
        }
    }
    return sqrt(sum);
}