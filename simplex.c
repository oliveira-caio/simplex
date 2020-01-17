#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

//essa função aloca um vetor de double com tamanho n para facilitar nossa vida durante o programa
void alocavetor(double **vetor, int n){

  int i;
  
  *vetor = (double*)malloc(n*sizeof(double));
  for(i=0; i<n; i++){
    (*vetor)[i] = 0;
  }

  return;
}

//essa função também serve para alocar um vetor de tamanho n, mas dessa vez vetores com entradas inteiras
void alocavetint(int **vetor, int n){

  int i;
  
  *vetor = (int*)malloc(n*sizeof(int));
  for(i=0; i<n; i++){
    (*vetor)[i] = 0;
  }
  
  return;
}

//essa função aloca uma matriz de double de tamanho m (linhas) por n (colunas) e preenche ela com zero
void alocamatriz(double ***matriz, int m, int n){

  int i, j;

  *matriz = (double**)malloc(m*sizeof(double*));
  
  for(i=0; i<m; i++){
    (*matriz)[i] = (double*)malloc(n*sizeof(double));
    for(j=0; j<n; j++){
      (*matriz)[i][j] = 0;
    }
  }

  return;
}

//essa função serve para desalocar uma matriz de double com m linhas
void desalocamatriz(double ***matriz, int m){

  int i;
  
  for(i=0; i<m; i++){
    free((*matriz)[i]);
    (*matriz)[i] = NULL;
  }
  
  free(*matriz);
  *matriz = NULL;

  return;
}

//função que desaloca vetores de inteiros
void desalocavetint(int **vetor){

  free(*vetor);
  *vetor = NULL;

  return;
}

//função que desaloca vetores de double
void desalocavetor(double **vetor){

  free(*vetor);
  *vetor = NULL;

  return;
}

//função que recebe uma matriz de double e inverte a matriz (nxn)
void invertematriz(double ***A, int n, int *erro){

  int i, j, k;
  double **aug2, x;

  alocamatriz(&aug2,n,2*n);
  if(aug2 == NULL){
    *erro = 2;
    return;
  }

  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      aug2[i][j] = (*A)[i][j];
    }
  }

  for(i=0; i<n; i++){
    for(j=n; j<2*n; j++){
      if(i == (j-n)){
	aug2[i][j] = 1;
      } else {
	aug2[i][j] = 0;
      }
    }
  }
  
  for(i=0; i<n; i++){
    if(aug2[i][i] == 0){
      for(j=i+1; j<n; j++){
	if(aug2[j][i] != 0){
	  for(k=0; k<2*n; k++){
	    x = aug2[i][k];
	    aug2[i][k] = aug2[j][k];
	    aug2[j][k] = x;
	  }
	}
      }
    }
  }
    
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      if(i != j){
	x = aug2[j][i]/aug2[i][i];
	for(k=0; k<2*n; k++){
	  aug2[j][k] -= x*aug2[i][k];
	}
      }
    }
  }
  
  for(i=0; i<n; i++){
    x = aug2[i][i];
    if(aug2[i][i] == 1){
      continue;
    } else {
      for(k=0; k<2*n; k++){
	aug2[i][k] = aug2[i][k]/x;
      }
    }
  }
  
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      (*A)[i][j] = aug2[i][j+n];
    }
  }
  
  desalocamatriz(&aug2,n);
  return;
  
}

//essa função recebe como parâmetros uma matriz M que será o resultado do produto das matrizes A e B. os inteiros a, b e c são, respectivamento, o número de linhas de A, o número de colunas de A/linhas de B e o número de colunas de B
void multmatriz(double ***M, double **A, double **B, int a, int b, int c){
  
  int i, j, k;
  
  for(i=0; i<a; i++){
    for(j=0; j<b; j++){
      (*M)[i][j] = 0;
      for(k=0; k<c; k++){
	(*M)[i][j] += A[i][k]*B[k][j];
      }
    }
  }

  return;
}

//função que vai de fato implementar o simplex. receberemos como parâmetro um vetor x, que será nosso vetor solução (passamos como **x para alterarmos o valor aqui na função); a nossa matriz A de restrições; o vetor b da igualdade das restrições (Ax = b); o vetor de custo c; uma solução inicial; os índices básicos; o número de linhas m da matriz A e o número de colunas n da matriz A; um contador para calcular o número de iterações; e por último uma variável de erro pra caso aconteça algo nas alocações
void simplex(double **x, double **A, double *b, double *c, double *solucao, int **base, int m, int n, int *cont, int *erro){
  
  int i, j, k = 0, l = -1, s, aux = 0, t;
  double cbarra, **cbt, **p, **B, **z, **u, thetaestrela, **Aj;
  
  *cont = *cont+1;
  if(*cont > 40000){
    *erro = 3;
    return;
  }
  alocamatriz(&cbt,1,m);
  alocamatriz(&Aj,m,1);
  alocamatriz(&p,1,m);
  alocamatriz(&B,m,m);
  alocamatriz(&z,1,1);
  u = NULL;

  if(cbt == NULL || Aj == NULL || p == NULL || B == NULL || z == NULL){
    *erro = 1;
    return;
  }
  
  for(i=0; i<m; i++){ 
    cbt[0][i] = c[(*base)[i]];
  }
  
  for(j=0; j<m; j++){
    for(i=0; i<m; i++){
      B[i][j] = A[i][(*base)[j]];
    }
  }

  for(i=0; i<m; i++){
    for(j=0; j<m; j++){
      printf("%lf ", B[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  invertematriz(&B,m,erro);
  if(*erro == 2){
    return;
  }
  multmatriz(&p,cbt,B,1,m,m);

  for(i=0; i<n; i++){
    for(j=0; j<m; j++){ 
      if(i == (*base)[j]){
	k++;
      }
    }
    if(k == 0){
      for(s=0; s<m; s++){
      	Aj[s][0] = A[s][i];
      }
      multmatriz(&z,p,Aj,1,1,m);
      cbarra = c[i] - **z;
      printf("custo: %.30lf\n", cbarra);
      if(cbarra < -pow(0.1,6)){
	alocamatriz(&u,m,1);
	if(u == NULL){
	  *erro = 1;
	  return;
	}
	multmatriz(&u,B,Aj,m,1,m);
	for(s=0; s<m; s++){
	  if(u[s][0] > pow(0.1,6)){
	    if(aux == 0){
	      thetaestrela = solucao[(*base)[s]]/u[s][0]; 
	      l = s;
	      aux++;
	    } else {
	      if(solucao[(*base)[s]]/u[s][0] < thetaestrela){
		thetaestrela = solucao[(*base)[s]]/u[s][0];
		l = s;
	      }
	    }
	  }
	}
	j = i;
	break;
      }
    }
    k=0;
  }

  if(u == NULL){
    for(i=0; i<n; i++){
      x[0][i] = solucao[i];
    }
    desalocamatriz(&cbt,1);
    desalocamatriz(&p,1);
    desalocamatriz(&B,m);
    desalocamatriz(&z,1);
    desalocamatriz(&Aj,m);
    return;
  } else {
    if(l == -1){
      desalocavetor(x);
      desalocamatriz(&cbt,1);
      desalocamatriz(&p,1);
      desalocamatriz(&B,m);
      desalocamatriz(&z,1);
      desalocamatriz(&u,m);
      desalocamatriz(&Aj,m);
      return;
    } else {
      for(s=0; s<m; s++){
	if(s != l){
	  solucao[(*base)[s]] = solucao[(*base)[s]] - thetaestrela*u[s][0];
	} else {
	  solucao[(*base)[s]] = 0;
	}
      }
      (*base)[l] = j;
      solucao[j] = thetaestrela;
    }
    desalocamatriz(&cbt,1);
    desalocamatriz(&p,1);
    desalocamatriz(&B,m);
    desalocamatriz(&z,1);
    desalocamatriz(&u,m);
    desalocamatriz(&Aj,m);
    simplex(x,A,b,c,solucao,base,m,n,cont,erro);
  }
}

int main(){

  int esc, esc1, m1, n1, m, n, i, j, k, *base, aux, tmp, cont = 0, erro;
  double **resposta80, **z, **P, **A, **B, **Aj, **Y, **D, **prod, *cy, *b, *c, *yini , *viavelcru, *viavel, *resposta, *aux2, aux3, soma = 0, tempo;
  clock_t t1;
  FILE *fp, *custo, *be, *game;
  
  srand(time(NULL));

  do{
  printf("Digite 1 se deseja encontrar uma estratégia ótima para o jogador-coluna em um jogo de soma-zero com dois jogadores.\nDigite 2 se deseja resolver um problema de otimização linear na forma padrão.\n");
  scanf("%d", &esc);
  } while(esc != 1 && esc != 2);

  if(esc == 1){
    do{
      printf("Digite 1 se deseja que a matriz de pagamento P seja lida de um arquivo.\nDigite 2 se deseja digitá-la.\n");
      scanf("%d", &esc1);
    } while(esc1 != 1 && esc1 != 2);
    
    if(esc1 == 1){
      fp = fopen("matriz.txt", "r");
      if(fp == NULL){
	printf("Arquivo com defeito.\n");
	return -1;
      }
      
      printf("Informe quantas linhas e quantas colunas tem a matriz de pagamento P: ");
      scanf("%d %d", &m1, &n1);
      
      alocamatriz(&P,m1,n1);
      if(P == NULL){
	printf("Erro de alocação.\n");
	return -1;
      } else {
	for(i=0; i<m1; i++){
	  for(j=0; j<n1; j++){
	    fscanf(fp,"%lf ",&(P[i][j]));
	  }
	}
	n = m1 + n1 + 2;
	m = m1 + 1;
      }
    }

    if(esc1 == 2){
      printf("Informe quantas linhas e quantas colunas tem a matriz de pagamento P: ");
      scanf("%d %d", &m1, &n1);
      
      alocamatriz(&P,m1,n1);
      if(P == NULL){
	printf("Erro de alocação.\n");
	return -1;
      } else {
	printf("Digite os elementos da matriz P:\n");
	for(i=0; i<m1; i++){
	  for(j=0; j<n1; j++){
	    printf("(%d,%d): ", i+1, j+1);
	    scanf("%lf", &(P[i][j]));
	  }
	}
	n = m1 + n1 + 2;
	m = m1 + 1;
      }
    }
  }
   
  if(esc == 2){ 
    printf("Considere um problema na forma padrão:\n\n     minimizar c^Tx\n     sujeito a Ax = b\n               x >= 0\n\n");
    do{
      printf("Informe quantas linhas e quantas colunas tem a matriz A: ");
      scanf("%d %d", &m, &n);
    } while((n < m) || (m <= 0) || (n <= 0));
  }

  alocamatriz(&z,1,1);
  alocavetor(&b,m);
  alocavetor(&c,n);
  alocamatriz(&A,m,n);
  alocamatriz(&B,m,m);
  alocamatriz(&Y,m,n+m);
  alocamatriz(&D,m,m);
  alocamatriz(&Aj,m,1);
  alocamatriz(&prod,m,1);
  alocavetor(&cy,m+n);
  alocavetor(&yini,m+n);
  alocavetor(&viavel,n);
  alocavetor(&viavelcru,n+m);
  alocavetor(&resposta,n);
  alocavetint(&base,m);

  if(b == NULL || c == NULL || A == NULL || B == NULL || Y == NULL || D == NULL || Aj == NULL || prod == NULL || cy == NULL || yini == NULL || viavel == NULL || viavelcru == NULL || resposta == NULL || base == NULL){
    printf("Erro de alocação.\n");
    return -1;
  } else {
    if(esc == 2){
      do{
      printf("Digite 1 se deseja que a matriz A de restrições seja lida de um arquivo.\nDigite 2 se deseja digitá-la.\n");
      scanf("%d", &esc1);
      } while(esc1 != 1 && esc1 != 2);
      
      if(esc1 == 1){
	fp = fopen("matriz.txt", "r");
	custo = fopen("vetc.txt", "r");
	be = fopen("vetb.txt","r");
	
	if(fp == NULL || custo == NULL || be == NULL){
	  printf("Arquivo com defeito.\n");
	  return -1;
	} else {
	  for(i=0; i<m; i++){
	    for(j=0; j<n; j++){
	      fscanf(fp,"%lf ", &(A[i][j]));
	    }
	  }
	  for(j=0; j<n; j++){
	    fscanf(custo,"%lf ", &(c[j]));
	  }
	  for(i=0; i<m; i++){
	    fscanf(be,"%lf ", &(b[i]));
	    if(b[i] < 0){
	      b[i] = -b[i];
	      for(j=0; j<n; j++){
		A[i][j] = -A[i][j];
	      }
	    }
	  }
	}
      }
	
      if(esc1 == 2){
	printf("Digite os elementos da matriz A:\n");
	for(i=0; i<m; i++){
	  for(j=0; j<n; j++){
	    printf("(%d,%d): ", i+1, j+1);
	    scanf("%lf", &(A[i][j]));
	  }
	}
	
	printf("Informe o vetor custo c: ");
	for(i=0; i<n; i++){
	  scanf("%lf", &c[i]);
	}
	
	printf("Informe o vetor b: ");
	for(i=0; i<m; i++){
	  scanf("%lf", &b[i]);
	  if(b[i] < 0){
	    b[i] = -b[i];
	    for(j=0; j<n; j++){
	      A[i][j] = -A[i][j];
	    }
	  }
	}
      }
    }
      
    if(esc == 1){
      for(j=0; j<n1; j++){
	A[0][j] = 1;
      }
      
      for(j=n1; j<n; j++){
	A[0][j] = 0;
      }
      
      for(i=1; i<m; i++){
	for(j=0; j<n1; j++){
	  A[i][j] = P[i-1][j];
	}
	
	A[i][n1] = 1;
	A[i][n1+1] = -1;
	
	for(j=n1+2; j<n; j++){
	  if(i == j-n1-1){
	    A[i][j] = 1;
	  } else {
	    A[i][j] = 0;
	  }
	}
      }
      
      c[n1] = -1;
      c[n1 + 1] = 1;
      b[0] = 1;
    }

    game = fopen("game.txt", "w");
    for(i=0; i<m; i++){
      for(j=0; j<n; j++){
	fprintf(game,"%.0lf ",A[i][j]);
      }
      fprintf(game,"0 0\n");
    }
    for(i=0; i<n; i++){
      for(j=0; j<n; j++){
    	if(i == j){
    	  fprintf(game,"%d ",1);
    	} else {
    	  fprintf(game,"%d ",0);
    	}
      }
      fprintf(game,"1 0\n");
    }
    
    for(i=0; i<m; i++){
      for(j=0; j<n; j++){
	Y[i][j] = A[i][j];
      }
    }
    
    for(i=0; i<m; i++){
      for(j=n; j<n+m; j++){
	if(i == j-n){
	  Y[i][j] = 1;
	} else {
	  Y[i][j] = 0;
	}
      }
    }
    
    for(i=0; i<n+m; i++){
      if(i < n){
	cy[i] = 0;
      } else {
	cy[i] = 1;
      }
    }
    
    for(i=0; i<n+m; i++){
      if(i < n){
	yini[i] = 0;
      } else {
	yini[i] = b[i-n];
      }
    }
    
    for(i=0; i<m; i++){
      base[i] = n+i;
    }
    
    t1 = clock();
    simplex(&viavelcru,Y,b,cy,yini,&base,m,n+m,&cont,&erro);
    if(erro == 3){
      printf("O número de iterações excedeu o limite.\n");
      return -1;
    } else if(erro == 1 || erro == 2){
      printf("Problema de alocação.\n");
      return -1;
    }
    t1 = clock() - t1;
    tempo = ((double)t1)/CLOCKS_PER_SEC;
    printf("\n\nTempo gasto no simplex do problema auxiliar: %lf segundos\n", tempo);
    printf("Número de iterações no problema auxiliar: %d\n", cont);
    cont = 0;
    
    for(i=0; i<n+m; i++){
      printf("%lf ", viavelcru[i]);
    }
    for(i=n; i<n+m; i++){
      soma += viavelcru[i];
    }
    
    if(soma > 0){
      printf("\n\nO problema é inviável.\n\n\n");
      return 0;
    }

    for(j=0; j<m; j++){
      for(k=0; k<m; k++){
	D[k][j] = Y[k][base[j]];
      }
    }
    invertematriz(&D,m,&erro);
    if(erro == 2){
      printf("Erro de alocação.\n");
      return -1;
    }

    for(i=0; i<m; i++){
      tmp = 0;
      if(base[i] >= n){
	for(j=0; j<n; j++){
	  for(k=0; k<m; k++){
	    Aj[k][0] = A[k][i];
	  }
	  multmatriz(&prod,D,Aj,m,1,m);
	  if(prod[i][0] != 0){
	    base[i] = j;
	    tmp = 1;
	    break;
	  } 
	}

	if(tmp == 0){
	  for(k=i; k<m-1; k++){
	    aux = base[k];
	    base[k] = base[k+1];
	    base[k+1] = aux;
	    aux2 = A[k];
	    A[k] = A[k+1];
	    A[k+1] = aux2;
	    aux3 = b[k];
	    b[k] = b[k+1];
	    b[k+1] = aux3;
	  }
	  m--;
	}
      }
    }
    
    for(i=0; i<n; i++){
      viavel[i] = viavelcru[i];
    }
    
    t1 = clock();
    simplex(&resposta,A,b,c,viavel,&base,m,n,&cont,&erro);
    if(erro == 3){
      printf("O número de iterações excedeu o limite.\n");
      return -1;
    } else if(erro == 1 || erro == 2) {
      printf("Erro de alocação.\n");
      return -1;
    }
    t1 = clock()-t1;
    tempo = (double)t1/CLOCKS_PER_SEC;
    printf("Tempo gasto no simplex do problema original foi: %lf seg\n", tempo);
    printf("Número de iterações no problema original: %d\n", cont);
    
    if(resposta == NULL){
      printf("\n\nO custo ótimo é menos infinito.\n");
      return 0;
    } else {
      if(esc == 2){
	printf("\n\n\nSolução ótima:\n");
	for(i=0; i<n; i++){
	  printf("x_%d = %lf\n", i+1, resposta[i]);
	}
	
	alocamatriz(&resposta80,n,1);
	if(resposta80 == NULL){
	  printf("Erro de alocação.\n");
	  return -1;
	}

	for(i=0;i<n;i++){
	  resposta80[i][0] = resposta[i];
	}
	multmatriz(&z,&c,resposta80,1,1,n);
	printf("Custo ótimo: %lf\n\n\n", **z);
	desalocamatriz(&resposta80,n);
      }
      
      if(esc == 1){
	printf("\n\nEstratégia ótima:\n");
	for(i=0;i<n1;i++){
	  printf("x_%d = %lf\n", i+1,resposta[i]);
	}
	printf("Ganho esperado: %lf\n\n\n", resposta[n1]-resposta[n1+1]);
      }
    }
  }

  desalocamatriz(&z,1);
  desalocavetor(&b);
  desalocavetor(&c);
  desalocamatriz(&A,m);
  desalocamatriz(&B,m);
  desalocamatriz(&Y,m);
  desalocamatriz(&D,m);
  desalocamatriz(&Aj,m);
  desalocamatriz(&prod,m);
  desalocavetor(&cy);
  desalocavetor(&yini);
  desalocavetor(&viavel);
  desalocavetor(&viavelcru);
  desalocavetor(&resposta);
  desalocavetint(&base);
  desalocamatriz(&P,m1);

  return 0;  
}

