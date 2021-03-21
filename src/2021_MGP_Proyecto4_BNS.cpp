#include <iostream>
using namespace std;

class Belief_Propagation
{
#define size_px 2
#define size_py 2
    float root[1][size_px]={0}; //Vector root

//Leafs, just for use
    float E[size_py][size_py]={0};
    float D[size_py][size_py]={0};
    float F[size_py][size_py]={0};

    float* L_E = 0;
    float* L_F = 0;
    float* L_D = 0;
    float* L_C;

    float* Pi_E=0;
    float* Pi_D_1=0;
    float* Pi_D=0;

    float evidence_vector[2] = {1,1};
    float default_evidence[2] = {1,1};

public:
    float* Posterior_E = 0;
    float* Posterior_C = 0;
    float* Posterior_D = 0;

    Belief_Propagation(int No_nodes, bool evidence, int evidence_case); //Constructors for 2 sizes of trees
    void Cases_generation(int evidence_case);
    float* Lambda_Operations(float Hyp_vector[], float *Prob);
    void Lambda_Operation(){
        L_E =new float[size_px];
        for(int i=0;i<size_px;i++)
            cout<<(L_E[i]=L_F[i]*L_D[i]);
    }
    float* Pi_Operation(float *Prob, float *Prob2);
    float* Post_Proba(float *Prob, float *Prob2);

};

int main() {
    Belief_Propagation b1(4, false, 2);


    return 0;
}
Belief_Propagation::Belief_Propagation(int No_nodes, bool evidence, int evidence_case) {
    if(evidence == false){
        evidence_vector[0] = 1;
        evidence_vector[1] = 0;
    }

    cout<<"Fill the CTPs for"<<No_nodes-1<<" and root"<<endl;
    cout<<"Root, CPT"<<endl;
    for(int i=0;i<size_px;i++)
        cin>>root[0][i];
    cout<<"Fill leafs"<<endl;
    if(No_nodes==4) {
        for (int i = 0; i < size_px; i++) {
            for (int j = 0; j < size_py; j++) {
                cin >> E[i][j];
            }
        } //1
        for (int i = 0; i < size_px; i++) {
            for (int j = 0; j < size_py; j++) {
                cin >> F[i][j];
            }
        } //2
        for (int i = 0; i < size_px; i++) {
            for (int j = 0; j < size_py; j++) {
                cin >> D[i][j];
            }
        } //3
    }
    else{
        cout<<"Sorry, We cant do this now"<<endl;
    }

    Cases_generation(evidence_case);

}
float* Belief_Propagation::Lambda_Operations(float Hyp_vector[], float *Prob) {
    float* Vector_Hyp = 0;
    Vector_Hyp =new float[size_px];
    float M[size_px][size_py]={0};
    for(int i=0;i<size_px;i++)
        for (int j = 0; j < size_py; j++)
            M[i][j]=*(Prob+i*size_px+j);

    for(int i=0;i<size_px;i++) {
        for (int j = 0; j < size_px; j++){
            Vector_Hyp[i] += Hyp_vector[j] * M[j][i];
        }
    }
    return Vector_Hyp;
    delete[] Vector_Hyp;

}
void Belief_Propagation::Cases_generation(int evidence_case) {
    if(evidence_case==2){
        cout<<"Evidence in root"<<endl;
        L_F = Lambda_Operations(evidence_vector, (float *) F);
        L_D = Lambda_Operations(default_evidence, (float *) D);

        Lambda_Operation();
        cout<<endl;
        L_C = Lambda_Operations(L_E,(float *)E);
        Pi_E = Pi_Operation((float*)E, (float*)root);

        Posterior_C = Post_Proba((float*)L_C,(float*)root);
        Posterior_E = Post_Proba((float*)L_E,(float*)Pi_E);

        Pi_D_1 =new float[size_px];
        cout<<endl;
        for(int i=0;i<size_px;i++)
            Pi_D_1[i]=Pi_E[i]*L_E[i];
        Pi_D = Pi_Operation((float*)D, (float*)Pi_D_1);

        Posterior_D = Post_Proba((float*)L_D,(float*)Pi_D);

    }

}
float* Belief_Propagation::Pi_Operation(float *Prob, float *Prob2) {
    float* Vector_Hyp = 0;
    Vector_Hyp =new float[size_px];

    float M[size_px][size_py]={0};
    float M1[1][size_py]={0};

    for(int i=0;i<size_px;i++)
        for (int j = 0; j < size_py; j++)
            M[i][j]=*(Prob+i*size_px+j);

    for(int i=0;i<1;i++)
        for (int j = 0; j < size_py; j++)
            M1[0][j]=*(Prob2+i*size_px+j);

    for(int i=0;i<size_px;i++) {
        for (int j = 0; j < size_px; j++){
            Vector_Hyp[i] += M1[0][j] * M[i][j];
        }
    }
    cout<<endl;
    return Vector_Hyp;
    delete[] Vector_Hyp;
}
float* Belief_Propagation::Post_Proba(float *Prob, float *Prob2) {
    float* Vector_Hyp = 0;
    Vector_Hyp =new float[size_px];
    float nom=0;

    float M[1][size_py]={0};
    float M1[1][size_py]={0};

    for(int i=0;i<1;i++)
        for (int j = 0; j < size_py; j++)
            M[0][j]=*(Prob+i*size_px+j);

    for(int i=0;i<1;i++)
        for (int j = 0; j < size_py; j++)
            M1[0][j]=*(Prob2+i*size_px+j);

    for(int i=0;i<size_px;i++) {
            Vector_Hyp[i] = M1[0][i] * M[0][i];
    }
    nom = 1/(Vector_Hyp[0] + Vector_Hyp[1]);
    for(int i=0;i<size_px;i++) {
        Vector_Hyp[i] = nom *Vector_Hyp[i];
    }
    cout<<endl;
    cout<<"Posterior_Probability"<<endl;
    for(int i=0;i<size_px;i++)
        cout<<Vector_Hyp[i]<<" ";
    cout<<endl;
    return Vector_Hyp;
    delete[] Vector_Hyp;
}