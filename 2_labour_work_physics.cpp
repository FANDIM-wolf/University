/*
 calculate main values in second labour work
 Result: variables for graph , which proves law of saving energy.

*/
#include <iostream>
#include <cmath>

using namespace std;




//1 formula
double calculate_translation_energy(double m2,double h, double t){
 //transform values
 m2 = m2 * pow(10,-3);
 h  = h * pow(10,-3);
 double translation_energy = (2 * m2 * (h*h))/t*t;

 return translation_energy ;

}
//2 formula
double calculate_translation_energy_Delta(double m2,double h, double t , double Dm2 , double Dh , double Dt){

    double translation_energy_Delta = (Dm2/m2)+(2*Dh/h)+(2*Dt/t);
    return translation_energy_Delta;

}
//3 formula
double calculate_translation_energy_D( double translation_energy_Delta , double translation_energy ){
     return translation_energy_Delta*translation_energy;
}
//4 formula
double calculate_rotational_energy(double h,double J0, double m , double R , double t , double d){
 //transform values
 h  = h * pow(10,-3);
 J0  = J0 * pow(10,-3);
 m =  m * pow(10,-3);
 R = R * pow(10,-3);
 d = d * pow(10,-3);

 double rotational_energy = 8*(h*h)*(J0+4*m*(R*R))/t*t*d*d;

 return rotational_energy ;

}
// 5 formula
double calculate_rotational_energy_Delta(double h,double J0, double m , double R , double t , double d , double Dh,double DJ0, double Dm , double DR , double Dt , double Dd){

    double rotational_energy_Delta_first_part_of_formula= (Dh/h)+(Dt/t)+(Dd/d); // values without transformation
    //transform values
    h  = h * pow(10,-3);
    J0  = J0 * pow(10,-3);
    m =  m * pow(10,-3);
    R = R * pow(10,-3);
    d = d * pow(10,-3);
    Dm = Dm * pow(10,-3);
    DR  = DR * pow(10,-3);
    double rotational_energy_Delta = 2 * rotational_energy_Delta_first_part_of_formula + ( (4*m*(R*R)*Dm ) + (8*m*R*DR) / J0+(4*m*(R*R)) );
    return rotational_energy_Delta;
}

//6 formula
double calculate_rotational_energy_D(double rotational_energy_Delta , double rotational_energy){

    return rotational_energy_Delta * rotational_energy;
}
//7 formula
double calculate_coefficient_of_friction(double m0 , double h){
    //transform values
    h  = h * pow(10,-3);
    m0 =  m0 * pow(10,-3);
    return  m0 * 10 *h;
}
// 8 formula
double calculate_coefficient_of_friction_Delta(double m0 , double h , double Dm0 , double Dh){

    return  (Dm0/m0)+(0.05/10)+(Dh/h);
}

// 9 formula
double calculate_coefficient_of_friction_D(double calculate_coefficient_of_friction , double calculate_coefficient_of_friction_Delta){

    return  calculate_coefficient_of_friction * calculate_coefficient_of_friction_Delta;
}

//10 formula
double calcualate_definite_translation_energy( double m2 , double h){
    //transform values
    h  = h * pow(10,-3);
    m2 =  m2 * pow(10,-3);
    return m2 * 10 * h;
}
//11 formula
double calcualate_definite_translation_energy_Delta( double m2 , double h , double Dm2 , double Dh){

    return  (Dm2/m2)+(0.05/10)+(Dh/h);
}
//12 formula
double calcualate_definite_translation_energy_D( double definite_translation_energy , double definite_translation_energy_D){

    return  definite_translation_energy * definite_translation_energy_D;
}
//13 formula
double calcualate_definite_translation_energy_with_comma(double translation_energy , double rotational_energy , double coefficient_of_friction){

 return translation_energy+rotational_energy+coefficient_of_friction;

}
//14 formula
double calcualate_definite_translation_energy_with_comma_Delta(double translation_energy_Delta , double rotational_energy_Delta , double coefficient_of_friction_Delta){
 return translation_energy_Delta + rotational_energy_Delta + coefficient_of_friction_Delta;
}
//15 formula
double calcualate_definite_translation_energy_with_comma_D(double definite_translation_energy_with_comma , double definite_translation_energy_with_comma_D){
 return definite_translation_energy_with_comma * definite_translation_energy_with_comma_D ;
}
int main()
{

    double  h , Dh ,d , Dd , m2 , Dm2 , R , DR , m , Dm , m0 , Dm0 , J0,DJ0 , t ,Dt ; // basic variables .
    cout << "2nd labour work in physics " << endl;
    cout << "v1.0 " << endl;
    cout <<"REMEBER! at moment when you are entering values with rational numbers , use '.' instead of  ',' and enter AVERAGE values .GOOD LUCK! "<<endl;
    //Enter basic variables
    cout<< "Height:";
    cin>>h;
    cout<< "Delta Height:";
    cin>>Dh ;
    cout<< "Diameter:";
    cin>>d ;
    cout<<"Delta Diameter:";
    cin>>Dd;
    cout<<"Mass of cylinder '(m2)' :";
    cin>>m2;
    cout<<"Delta Mass of cylinder '(m2)' :";
    cin>>Dm2;
    cout<<"Radius :";
    cin>>R;
    cout<<"Delta Radius:";
    cin>>DR;
    cout<<"Mass of small cylinder '(there are four in total)': ";
    cin>>m;
    cout<<"Delta Mass of small cylinder '(there are four in total)': ";
    cin>>Dm;
    cout<<"Enter mass 0 '(m0)'";
    cin>>m0;
    cout<<"Enter Delta mass 0 '(m0):";
    cin>>Dm0;
    cout<<"Enter moment of Inertion J0 : ";
    cin>>J0;
    cout<<"Enter Delta moment of Inertion J0 :";
    cin>>DJ0;
    cout<<"Enter time '(t)':";
    cin>>t;
    cout<<"Enter Delta time:";
    cin>>Dt;
    //first Test :: print some values
    cout<<"Small test , printing of some values"<<endl;
    cout<<h<<Dh<<d<<Dd<<m2<<Dm2<<"Momemt of Inertion:"<< J0 <<endl;

    //main part
    double translation_energy = calculate_translation_energy(m2,h,t);
    cout<<"Translation energy:"<<translation_energy<<endl;
    double translation_energy_Delta = calculate_translation_energy_Delta(m2,h,t,Dm2,Dh,Dt);
    cout<<"Translation energy Delta:"<<translation_energy_Delta<<endl;
    double translation_energy_D = calculate_translation_energy_D(translation_energy_Delta,translation_energy);
    cout<< "Translation energy D: "<<translation_energy_D<<endl;
    double rotational_energy = calculate_rotational_energy(h,J0,m,R,t ,d);
    cout<<"Rotational energy:"<<rotational_energy<<endl;
    double rotational_energy_Delta = calculate_rotational_energy_Delta(h,J0,m,R,t,d,Dh,DJ0,Dm,DR,Dt,Dd);
    cout<<"Rotational energy Delta:"<<rotational_energy_Delta<<endl;
    double rotational_energy_D = calculate_rotational_energy_D(rotational_energy , rotational_energy_Delta);
    cout<< "Rotational energy D :"<<rotational_energy_D<<endl;
    double coefficient_of_friction = calculate_coefficient_of_friction(m0,h);
    cout<<"Coefficient of friction :"<<coefficient_of_friction<<endl;
    double coefficient_of_friction_Delta = calculate_coefficient_of_friction_Delta(m0,h,Dm0,Dh);
    cout<<"Coefficient of friction Delta :"<<coefficient_of_friction_Delta<<endl;
    double coefficient_of_friction_D = calculate_coefficient_of_friction_D(coefficient_of_friction ,coefficient_of_friction_Delta);
    cout<<"Coefficient of friction D:" << coefficient_of_friction_D <<endl;
    double definite_translation_energy = calcualate_definite_translation_energy(m2,h);
    cout<<"Definite translation energy:"<<definite_translation_energy;
    double definite_translation_energy_Delta = calcualate_definite_translation_energy_Delta(m2,h,Dm2,Dh);
    cout<<"Definite translation energy Delta: "<<definite_translation_energy_Delta<<endl;
    double definite_translation_energy_D = calcualate_definite_translation_energy_D(definite_translation_energy ,definite_translation_energy_Delta);
    cout<< "Definite translation energy D "<<definite_translation_energy_D<<endl;


    double definite_translation_energy_with_comma = calcualate_definite_translation_energy_with_comma(translation_energy,rotational_energy,coefficient_of_friction);
    cout<< "Definite translation energy with comma:"<<definite_translation_energy_with_comma<<endl;

    double definite_translation_energy_with_comma_Delta = calcualate_definite_translation_energy_with_comma_Delta(translation_energy_D , rotational_energy_D , coefficient_of_friction_D);
    cout<< "Definite translation energy with comma Delta:"<<definite_translation_energy_Delta<<endl;

    double definite_translation_energy_with_comma_D = calcualate_definite_translation_energy_with_comma_D(definite_translation_energy_with_comma  ,definite_translation_energy_with_comma_Delta);
    cout<< "Definite translation energy with comma D"<<definite_translation_energy_with_comma_D<<endl;


    return 0;
}
