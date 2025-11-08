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
    free(Av);
    free(z);
}

void topkSingular(double **A, int m, int n, int k, double **U, double **V, double *sing){//computing top k singular values
    if(k <=0) return;
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

    for(int i =0;i<k; i++){// k singular values joiye che
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

        if(sigma >1e-12){// u (A*v) /sigma;
            vec_devide(Av, m, sigma);
        }else{
            for(int j =0;j<m; j++){ 
                Av[j] =0.0;
            }
        }
        
        //u & v ne U & V ma nakhva
        for(int j =0;j<n; j++){
            V[i][j] = v[j];// V(k n) - This is V_kt
        }
        for(int j =0;j<m; j++){
            U[j][i] = Av[j];// U(m k)
        }

        if(sigma >1e-12){
            for(int j =0; j<m; j++){//Acopy ne deflate karo
                for(int l =0; l<n; l++){// k thi confusion no thay
                    A_copy[j][l] = A_copy[j][l] - sing[i]*U[j][i]*V[i][l];//deflate
                }
            }
        }
    }

    // All printf and free calls removed. Main will handle this.
    free_mat(A_copy, m);
    free(v);
    free(Av);
}

int read_image_rgb(const char *filename, int *width,int *height, double **R, double **G, double **B){
    int n_channels;
    unsigned char *img_data = stbi_load(filename, width,height, &n_channels, 3);//n_channel ni value pan 3 che(RGB)
    
    if(img_data == NULL){
        fprintf(stderr, "Error: Failed to load image %s. Reason: %s\n", filename, stbi_failure_reason());
        return 0;
    }
    printf("Successfully read %s.\n", filename);

    for(int y =0;y<*height; y++){//Copy kare data unsigned char buffer thi double matrices
        for(int x =0;x<*width; x++){
            unsigned char *pixel = img_data+(y*(*width) + x)*3;//pointer ne pixel na starting(y, x) ma leva
            R[y][x] =(double)pixel[0];//Red
            G[y][x] =(double)pixel[1];//Green
            B[y][x] =(double)pixel[2];//Blue
        }
    }
    stbi_image_free(img_data);//Loaded data free kare ,stbi_image thi natu thatu
    return 1;
}

void normalize_channel(double **matrix, int width,int height, unsigned char *buffer, int offset, int stride){
    double min_val =matrix[0][0];
    double max_val =matrix[0][0];
    for(int y =0;y<height;y++){
        for(int x =0;x<width; x++){
            if(matrix[y][x] < min_val){
                min_val = matrix[y][x];
            }
            if(matrix[y][x] > max_val){
                max_val = matrix[y][x];
            }
        }
    }

    double range = max_val - min_val;
    if(range < 1e-9){//Handle flat matrix
        range =1.0; 
    }

    for(int y =0;y<height; y++){
        for(int x =0;x<width; x++){
            double norm_val = (matrix[y][x] - min_val)/range;
            if(norm_val < 0.0){// Clamp values
                norm_val = 0.0;
            }
            if(norm_val > 1.0){
                norm_val = 1.0;
            }
            buffer[(y * width + x) * stride + offset] = (unsigned char)(norm_val * 255.0);//buffer ma sachi jagya a nakhvu
        }
    }
}

void write_image_rgb(const char *filename, int width,int height, int quality, double **R,double **G,double **B){
    //buffer allocate karyu aakhi image mate (width * height * 3 channels);
    unsigned char *img_data =(unsigned char*)malloc(width*height*3*sizeof(unsigned char));
    if(img_data == NULL){
        fprintf(stderr,"Error: Malloc failed for image buffer.\n");
        return;
    }

    //Normalize each channel and place it in the interleaved buffer
    normalize_channel(R, width, height, img_data, 0, 3); // Offset 0, Stride 3 [red]normalizing Red channel
    normalize_channel(G, width, height, img_data, 1, 3); // Offset 1, Stride 3 [green]
    normalize_channel(B, width, height, img_data, 2, 3); // Offset 2, Stride 3 [blue]

    // stbi_write_jpg(filename, width, height, channels, data, quality)
    if(stbi_write_jpg(filename, width,height, 3, img_data, quality) == 0){
        fprintf(stderr,"Error: Failed to write image %s\n",filename);
    }else{
        printf("Successfully wrote matrix to %s.\n",filename);
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

int main(){
    int m,n;//(height, width)
    
    const char* filename = "inputbhai.jpg";// jpg laine chalu chu

    int n_channels_check;//matrices allocate kariye R,G,B channel mate; ena mate m n known hova joiye creat_mat mate to pela size lai laiye
    if(!stbi_info(filename, &n, &m, &n_channels_check)){
        filename = "inputbhai.png";//jo jpg nathi to png try kariye
        printf("Couldn't find .jpg, trying .png\n");
        if(!stbi_info(filename, &n, &m, &n_channels_check)){// == NULL karva thi comparision bw ptr and int aavi jay che
            fprintf(stderr,"Failed to read info for input image (tried .jpg and .png).\n");
            return 1;
        }
    }

    double **R =creat_mat(m, n);// have m n aavi gya
    double **G =creat_mat(m, n);
    double **B =creat_mat(m, n);

    if(!read_image_rgb(filename, &n, &m, R, G, B)){//have image ne matrix ma read karo
        fprintf(stderr,"Failed to read input image.\n");
        return 1;
    }

    int k;// k input lyo
    printf("Enter the number of singular values (k): ");
    scanf("%d",&k);
    int k_max = (n < m) ? n : m;//k max rank thi toh nano j hovo joiye max rnak = min(dimension)
    if(k>k_max){
        printf("Changing k (%d) to max possible (%d).\n", k, k_max);
        k = k_max;
    }
    if(k <=0){
        printf("Error: rank(k) must be positive.\n");
        return 1;
    }

    double **R_k =creat_mat(m, n);//reconstructed matrix mate allocation
    double **G_k =creat_mat(m, n);
    double **B_k =creat_mat(m, n);

    //red badha mate svd vaprvu pade
    double **U_r =creat_mat(m, k);
    double **V_r =creat_mat(k, n);
    double *sing_r =malloc(k*sizeof(double));
    
    printf("Calculating top %d singular values (RED)\n", k);
    topkSingular(R, m,n, k, U_r,V_r, sing_r);
    
    for(int i =0;i<m; i++){//topk U_r V_r singular value aape ene matrix ma store karva
        for(int j =0;j<n; j++){
            double sum =0.0;
            for(int l =0;l<k; l++){
                sum += U_r[i][l]*sing_r[l]*V_r[l][j];
            }
            R_k[i][j] =sum;
        }
    }

    free_mat(U_r, m); free_mat(V_r, k); free(sing_r);
    printf("Red Done.\n");

    //green
    double **U_g =creat_mat(m, k);
    double **V_g =creat_mat(k, n);
    double *sing_g =malloc(k*sizeof(double));
    
    printf("Calculating top %d singular values (GREEN)...\n", k);
    topkSingular(G, m,n, k, U_g, V_g, sing_g);
    
    for(int i =0;i< m; i++){
        for(int j =0;j<n; j++){
            double sum =0.0;
            for(int l =0;l<k; l++){
                sum += U_g[i][l]*sing_g[l]*V_g[l][j];
            }
            G_k[i][j] =sum;
        }
    }

    free_mat(U_g, m); free_mat(V_g, k); free(sing_g);
    printf("Green Done.\n");

    //blue
    double **U_b =creat_mat(m, k);
    double **V_b =creat_mat(k, n);
    double *sing_b =malloc(k*sizeof(double));
    
    printf("Calculating top %d singular values (BLUE)...\n", k);
    topkSingular(B, m,n, k, U_b, V_b, sing_b);
    
    for(int i =0;i<m; i++){
        for(int j =0;j<n; j++){
            double sum =0.0;
            for(int l =0;l<k; l++){
                sum += U_b[i][l]*sing_b[l]*V_b[l][j];
            }
            B_k[i][j] =sum;
        }
    }

    free_mat(U_b, m); free_mat(V_b, k); free(sing_b);
    printf("Blue Done.\n");

    printf("\nImage done\n");// traney ne RGB vala write kari gyo
    write_image_rgb("output_test.jpg", n,m, 90, R_k, G_k, B_k);//90 quality
    
    printf("freeing memory.\n");
    free_mat(R, m); free_mat(G, m); free_mat(B, m);
    free_mat(R_k, m); free_mat(G_k, m); free_mat(B_k, m);

    printf("Done.\n");
    return 0;
}