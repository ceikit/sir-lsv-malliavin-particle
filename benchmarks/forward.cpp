// Conservative forward generator on a uniform (log S, x_d, x_f) grid.
// Correlated diffusion uses diagonal jumps; central axial drift is retained
// wherever positive and otherwise replaced by the minimal positive pair with
// the same drift. The latter adds numerical diffusion in remote drift regions.
#include <algorithm>
#include <cmath>
#include <cstddef>

extern "C" double rhs(const double* y, double* out, int nf, int nx, int nu, int nv,
                       double hx, double hu, double hv, const double* sig,
                       const double* slope, const double* par, double t) {
    const std::size_t N = std::size_t(nx)*nu*nv;
    std::fill(out, out+nf*N, 0.0);
    const double kd=par[0], sd=par[1], kf=par[2], sf=par[3];
    const double beta=sd*par[4], gamma=sf*par[5], c=sd*sf*par[6];
    const auto B=[](double k, double t){return k==0 ? t : -std::expm1(-k*t)/k;};
    const double cd=.5*sd*sd*std::pow(B(kd,t),2), cf=.5*sf*sf*std::pow(B(kf,t),2);
    const double vd=sd*sd*B(2*kd,t), vf=sf*sf*B(2*kf,t);
    const double h[3]={hx,hu,hv};
    const std::ptrdiff_t stride[3]={nu*nv,nv,1};
    const int sizes[3]={nx,nu,nv};
    double maxrate=0;
    for(int i=0;i<nx;++i) {
        const double s=sig[i], g=slope[i];
        const double q[3]={s*s,sd*sd,sf*sf};
        const double cross[3]={beta*s,gamma*s,c};
        double axial[3]={q[0],q[1],q[2]};
        const int aa[3]={0,0,1}, bb[3]={1,2,2};
        for(int pair=0;pair<3;++pair){
            int a=aa[pair], b=bb[pair];
            axial[a]-=std::abs(cross[pair])*h[a]/h[b];
            axial[b]-=std::abs(cross[pair])*h[b]/h[a];
        }
        if(*std::min_element(axial,axial+3)<-1e-14) return -1.0;
        for(int j=0;j<nu;++j) for(int k=0;k<nv;++k) {
            const int index[3]={i,j,k};
            const std::size_t n=(std::size_t(i)*nu+j)*nv+k;
            const double u=(j-nu/2)*hu, v=(k-nv/2)*hv;
            const double rd=par[7]+cd+u;
            const double drift[3]={rd-par[8]-cf-v-.5*s*s,-kd*u,-kf*v-gamma*s};
            double total=0;
            auto jump=[&](std::ptrdiff_t offset,double rate){
                if(rate==0)return;
                total+=rate;
                for(int f=0;f<nf;++f) {
                    double flux=rate*y[f*N+n];
                    out[f*N+n+offset]+=flux;
                    out[f*N+n]-=flux;
                }
            };
            for(int a=0;a<3;++a){
                const double adv=drift[a]/h[a];
                const double diffusion=std::max(axial[a]/(2*h[a]*h[a]),.5*std::abs(adv));
                if(index[a]+1<sizes[a])jump(stride[a],diffusion+.5*adv);
                if(index[a]>0)jump(-stride[a],diffusion-.5*adv);
            }
            for(int pair=0;pair<3;++pair){
                const int a=aa[pair], b=bb[pair], sign=cross[pair]>=0?1:-1;
                const double rate=std::abs(cross[pair])/(2*h[a]*h[b]);
                for(int d : {-1,1}){
                    if(index[a]+d>=0 && index[a]+d<sizes[a] &&
                       index[b]+d*sign>=0 && index[b]+d*sign<sizes[b])
                        jump(d*(stride[a]+sign*stride[b]),rate);
                }
            }
            maxrate=std::max(maxrate,total+std::max(rd,0.0));
            for(int f=0;f<nf;++f) out[f*N+n]-=rd*y[f*N+n];
            if(nf==5){
                const double p=y[n], m1=y[N+n],m2=y[2*N+n],m3=y[3*N+n],m4=y[4*N+n];
                out[N+n]+=-(kd+s*g)*m1-m2+(vd+beta*s)*p;
                out[2*N+n]+=-(kd+kf)*m2-gamma*g*m1+c*p;
                out[3*N+n]+=-kf*m3+m4+(-vf+gamma*s)*p;
                out[4*N+n]+=-(kd+kf)*m4+beta*g*m3+c*p;
                // Exactly the shared-noise flux -div(g h m), h=(s,beta,gamma).
                const double gh[3]={g*s,g*beta,g*gamma};
                for(int f : {1,3}) for(int a=0;a<3;++a) for(int d : {-1,1}){
                    if(index[a]+d>=0 && index[a]+d<sizes[a]){
                        const double flux=d*gh[a]*y[f*N+n]/(2*h[a]);
                        out[f*N+n+d*stride[a]]+=flux;
                        out[f*N+n]-=flux;
                    }
                }
            }
        }
    }
    return maxrate;
}
