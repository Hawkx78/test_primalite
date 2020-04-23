#include "stdio.h"
#include "stdlib.h"
#include "gmp.h"


void pgcd(mpz_t n, mpz_t d, mpz_t b){

	mpz_t a;
	mpz_t r;
	mpz_t r1;
	
	mpz_init(a);
	mpz_init(r);
	mpz_init(r1);
	
	mpz_set(a,n);
	mpz_set(b,d);

	mpz_mod(r, a, b);
	if(mpz_cmp_ui( r, 0) == 0){ mpz_set(b,d);}
	else {
		mpz_mod(r1, a, b);
		pgcd(b, r1, b);
	} 
}


/*int JacobiF(mpz_t num, mpz_t den,int *sign){
	
	mpz_t b;
	mpz_t p;
	mpz_t numt;
	mpz_t dent;
	mpz_t tmp;
	
	mpz_init(b);
	mpz_init(p);
	mpz_init(numt);
	mpz_init(dent);
	mpz_init(tmp);
	
	//gmp_printf(" num = %Zd", num);
	//gmp_printf(" den = %Zd\n", den);
	if (*sign > 0)
	  gmp_printf("%Zd/%Zd=\n", num, den);
	else
	  gmp_printf("-(%Zd/%Zd)=\n", num, den);
	
	mpz_mod(num, num, den);
	
	gmp_printf(" num = %Zd", num);
	gmp_printf(" den = %Zd\n", den);
	
	mpz_mod_ui(tmp ,num, 2);
	mpz_set(numt,num);
	gmp_printf(" num = %Zd", numt);
	gmp_printf(" tmp = %Zd\n", tmp);
	
	for(mpz_set_str(p, "1", 10); mpz_cmp_ui(tmp, 0) == 0 ; mpz_add_ui(p , p, 1)){
		if(mpz_cmp_ui(numt, 2) > 0){
	
			mpz_cdiv_q_ui(numt, numt, 2);
			mpz_mod_ui(tmp ,numt, 2);
			gmp_printf(" num1 = %Zd", numt);
			gmp_printf(" tmp1 = %Zd\n", tmp);
			mpz_set(num, numt);
		}else{
	
			mpz_set_str(tmp, "9", 10);
			mpz_set(num, numt);
			gmp_printf(" pui = %Zd\n", p);
		}
	}
	
	gmp_printf(" numV = %Zd", num);
	gmp_printf(" denV = %Zd", den);
	if(mpz_cmp_ui(num, 1) == 0){
		return 1;
	}
	if(mpz_cmp_ui(num, 2) == 0){
		mpz_mod_ui(dent, den, 8);
		gmp_printf(" dent = %Zd\n", dent);
		if(mpz_cmp_ui( dent, 7) == 0 || mpz_cmp_ui( dent, 1) == 0) return 1;
		else if(mpz_cmp_ui( dent, 3) == 0 || mpz_cmp_ui( dent, 5) == 0) 
		  { return -1; }
	}
	
	pgcd(num, den, b);
	gmp_printf(" b = %Zd\n", b);
	if(mpz_cmp_ui( b, 1) == 0){
		
		mpz_set(numt,num);
		mpz_set(dent,den);
		mpz_set(tmp,den);
		mpz_mod_ui(numt, num, 4);
		mpz_mod_ui(dent, den, 4);
		gmp_printf(" numt = %Zd", numt);
		gmp_printf(" dent = %Zd\n", dent);
		if(mpz_cmp_ui( numt, 3) == 0 && mpz_cmp_ui( dent, 3) == 0){
			
			mpz_set(den,num);
			mpz_set(num,tmp);
			 *sign = *sign * -1;
			
			gmp_printf(" num = %Zd", num);
			gmp_printf(" den = %Zd\n", den);
			return JacobiF(num, den, sign);
		}
		else if(mpz_cmp_ui( numt, 1) == 0 || mpz_cmp_ui( dent, 1) == 0){
			
			mpz_set(den,num);
			gmp_printf(" tmp = %Zd", tmp);
			mpz_set(num,tmp);
			gmp_printf(" num = %Zd", num);
			gmp_printf(" den = %Zd\n", den);
			return JacobiF(num, den, sign);
		}
	}else {
		return 0;
	}
	return 99999;
}*/

/*int Jacobi(mpz_t num, mpz_t den) {
	int sign = +1;
	int res = JacobiF(num, den, &sign);
	printf("%d", res * sign);
	return res * sign;
}*/


int Jacobi(mpz_t num, mpz_t den, mpz_t signe){
	
	mpz_t tmp;
	mpz_t numtmp;
	mpz_t dentmp;
	mpz_t p;
	mpz_t Respgcd;
	
	mpz_init(p);
	mpz_init(numtmp);
	mpz_init(dentmp);
	mpz_init(tmp);
	mpz_init(Respgcd);
	
	//propriété 1
	mpz_mod(num, num, den);
	gmp_printf("Propriété 1 Num = %Zd, Den = %Zd, Signe = %Zd\n", num, den, signe);
	
	//Propriété 3
	mpz_mod_ui(tmp, num, 2);
	for(mpz_set_str(p, "0", 10); mpz_cmp_ui(tmp, 0) == 0 ; mpz_add_ui(p , p, 1)){
		if(mpz_cmp_ui(num, 1) > 0){
			mpz_cdiv_q_ui(numtmp, num, 2);
			mpz_mod_ui(tmp ,numtmp, 2);
			mpz_set(num, numtmp);
			gmp_printf("Propriete 3 Num = %Zd Den = %Zd", num, den);
			gmp_printf(" Signe = %Zd ", signe);
			
			//Propriété 5
			{
				gmp_printf(" puissance = %Zd ", p);
				mpz_mod_ui(dentmp, den, 8);
				gmp_printf(" Mod = %Zd ", dentmp);
				if(mpz_cmp_ui( dentmp, 7) == 0 || mpz_cmp_ui( dentmp, 1) == 0){ // si n=1 mod 8 ou n =7 mod 8
					//mpz_neg(signe, signe);
					gmp_printf("Propriete 5.1 Signe = %Zd\n",signe);; //return 1;
				}
				else if(mpz_cmp_ui( dentmp, 3) == 0 || mpz_cmp_ui( dentmp, 5) == 0){ // si  n=3 mod 8 ou n =5 mod 8
						mpz_neg(signe, signe);
						gmp_printf("Propriete 5.2 Signe = %Zd\n",signe);
				}
			}
			
		}else if(mpz_cmp_ui(num, 2) == 0){
			
			mpz_mod_ui(p, p, 2); // Pour faire le calcul avec la puissance
			if(mpz_cmp_ui(p, 0) == 0){
				if(mpz_cmp_ui(signe, -1)){
					mpz_neg(signe, signe);
				}
				gmp_printf("Propriete 5.3 Signe = %Zd\n",signe); gmp_printf("Propriete 5.3.2 Num = %Zd Den = %Zd", num, den);return 1;
			}else{
				if(mpz_cmp_ui(signe, 1)){
					mpz_neg(signe, signe);
				}
				gmp_printf("Propriete 5.4 Signe = %Zd\n",signe);gmp_printf("Propriete 5.4.2 Num = %Zd Den = %Zd", num, den); return 1;
			}
		}else{
			mpz_set_str(tmp, "9", 10);
			mpz_set(num, numtmp);
			gmp_printf("Propriete 3 Num = %Zd, Den = %Zd", num, den);
		}
	}

	//Propriété 4
	if(mpz_cmp_ui( num, 1) == 0){/*
		if(mpz_cmp_ui(signe, -1)){
			gmp_printf("Propriete 4 Res = %Zd\n", signe); return -1;
		}else{*/
			gmp_printf("Propriete 4 Res = %Zd\n", signe); return 1;
		/*}*/
	}
	
	//Propriete 2
	pgcd(num, den, Respgcd); 
	gmp_printf("Propriete 2 Num = %Zd, Den = %Zd, pgcd = %Zd\n", num, den, Respgcd); 
	if(mpz_cmp_ui( Respgcd, 1) != 0){
		gmp_printf("Propriete 2 Res = 0\n"); return 0;
	} else {
		
		//Propriete 6
		//On a m/n
		mpz_set(numtmp,num);
		mpz_set(dentmp,den);
		
		mpz_mod_ui(numtmp, num, 4);
		mpz_mod_ui(dentmp, den, 4);
		if(mpz_cmp_ui( numtmp, 3) == 0 && mpz_cmp_ui( dentmp, 3) == 0){// si n=3 mod 4 et n =3 mod 4
			//Ici on doit faire -(n/m) donc on echange 
			mpz_set(tmp,num);
			mpz_set(num,den);
			mpz_set(den,tmp);
			mpz_neg(signe, signe);
			gmp_printf("Propriete 6 Num = %Zd, Den = %Zd, Signe = %Zd\n", num, den, signe); 
			return Jacobi(num, den, signe); // Appel récursif pour revenir à l'étape 1
		}
		else if(mpz_cmp_ui( numtmp, 1) == 0 || mpz_cmp_ui( dentmp, 1) == 0){// si n=1 mod 4 ou m =1 mod 4
			//Ici on doit faire n/m
			mpz_set(tmp,num);
			mpz_set(num,den);
			mpz_set(den,tmp);
			gmp_printf("Propriete 6 Num = %Zd, Den = %Zd\n", num, den); 
			return Jacobi(num, den, signe); //Appel récursif pour revenir a l'étape 1
		}
	}
	return 9999; //Return si erreur
}

void square_and_multiply( mpz_t a, mpz_t n, mpz_t H){
	
	mpz_t r;
	mpz_init(r);
	mpz_set(r,a); // r = a

	gmp_printf("\n a = %Zd", a);
	gmp_printf("\n r = %Zd", r);
	
	
	for (int i = 0; mpz_get_str(NULL, 2, H)[i]-48 != -48 ; i++)
	{
		printf("\n%d", mpz_get_str(NULL, 2, H)[i]-48);
	}
	printf("\n");
	
	for (int i = 1; mpz_get_str(NULL, 2, H)[i]-48 != -48; i++){
			
			mpz_mul( r, r, r);
			mpz_mod( r, r, n);  //r = (r * r) % n;
			gmp_printf("\n a2 = %Zd", a);
			gmp_printf("\n r2 = %Zd\n", r);
		if ((mpz_get_str(NULL, 2, H)[i]-48) == 1){
			
			mpz_mul( r, r, a); 
			mpz_mod( r, r, n);//r = r * a % n
			gmp_printf("\n a3 = %Zd", a);
			gmp_printf("\n r3 = %Zd\n", r);
		}
	}
	mpz_set(a,r);
	gmp_printf("\nres = %Zd\n\n", r);
}

int Solovay_Strassen(mpz_t n, mpz_t k){
	
	int r = 999;
	
	mpz_t a;
	mpz_t atmp;
	mpz_t ntmp;
	mpz_t rtmp;
	mpz_t rt;
	mpz_t s;
	mpz_t signe;
	
	mpz_init(a);	
	mpz_init(atmp);	
	mpz_init(ntmp);	
	mpz_init(rtmp);
	mpz_init(rt);
	mpz_init(s);
	mpz_init(signe);
	
	mpz_set_str(signe,"1",10);
	
	gmp_randstate_t state;
	gmp_randinit_mt(state); 
	
	if (mpz_cmp_ui( n, 2) < 0) return 0; 
	mpz_mod_ui(ntmp, n, 2);
    if (mpz_cmp_ui( n, 2) != 0 && mpz_cmp_ui( ntmp, 0) == 0) return 0; 
	
	for (int i = 0; mpz_cmp_ui( k, i) > 0; i++){
		mpz_set_str(signe,"1",10);
		mpz_sub_ui(ntmp, n, 1);
		mpz_urandomm( a, state, n);
		while(mpz_cmp_ui(a, 0) == 0 || mpz_cmp_ui(a, 1) == 0 ) {mpz_urandomm( a, state, n);gmp_printf("RAND = %Zd,", a); printf("i = %d\n",i);}
			
		gmp_printf("RAND = %Zd,", a); printf("i = %d\n",i);
		mpz_set(atmp, a);
		mpz_set(ntmp, n);
		r = Jacobi(	atmp, ntmp, signe);
		gmp_printf("signe =%Zd\n",signe);
		mpz_mul_ui(rt, signe, r);
		gmp_printf("rt =%Zd\n",rt);
		
		mpz_sub_ui( s, n, 1);
		mpz_cdiv_q_ui( s, s, 2); 
		
		square_and_multiply( a, n, s);
		mpz_mod(rtmp, rt, n);
		gmp_printf("A = %Zd, rtmp =%Zd, n = %Zd\n", a, rtmp, n);
		if( r == 0 || mpz_cmp( a, rtmp) != 0){
			return 0; //composé
		}
	}
	return 1; //premier
}

int main(int argc, char *argv[]){
	
	mpz_t k;
	mpz_t n;
	mpz_t signe;
	
	mpz_init(n);
	mpz_init(k);
	mpz_init(signe);
	
	
	int i = 100;
	int res = -1;
	
	mpz_set_str(n, argv[1], 10);
	mpz_set_str(k, argv[2], 10);
	mpz_set_str(signe, "1", 10);

	
	/*square_and_multiply(a,b,c);*/
	
	//~ i = Jacobi(n, k, signe);
	//~ printf("i=%d, signe=", i);
	//~ gmp_printf("%Zd\n", signe);
	
	res = Solovay_Strassen(n,k);
	if(res == 0) printf("composé\n");
	else if(res == 1) printf("\npremier\n");
	
	/*printf("Entrer votre nombre : ");
	gmp_scanf("%Zd\n", n);
	//scanf("%d\n", &k);
	//mpz_set_str(n, "296217028895100596785704895132502333533045130165114577343448265793942220615707859996559248687671450426890276572374230621654480300851501555796620774175002000585574362622003314933153119011598058864406732063707925937471648439951799503525544290538150979346898922897057835980982976679408923392726902293797105950859", 10);
	gmp_printf("nombre = %Zd\n", &n);
	//printf("repet = %d\n", k);*/
	
	mpz_clear(n);
	mpz_clear(k);
	
}
