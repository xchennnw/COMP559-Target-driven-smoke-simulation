package comp559.fluid;

import java.util.LinkedList;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2f;
import javax.vecmath.Tuple2f;

import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;


/**
 * Eulerian fluid simulation class. 
 * 
 * This follows closely the solver documented in J. Stam GDC 2003, specifically
 * this uses a non staggered grid, with an extra set of grid cells to help enforce
 * boundary conditions
 * 
 * @author kry
 */
public class Fluid {
    
    final static int DIM = 2;
    
    /** 
     * Velocity
     * U[0][IX(0,0)] is the x velocity in the (0,0) grid location.
     * */ 
    public float[][] U0;
    
    /** temporary velocity variable */
    private float[][] U1;          
       
   
    private IntParameter Nval = new IntParameter( "grid size", 100, 4, 22500);
    
    /** Number of grid cells (not counting extra boundary cells */
    public int N = 150;
    
    /** Dimension of each grid cell */
    public float dx = 1;
    public float dt = 0.1f;
    
    /** time elapsed in the fluid simulation */
    public double elapsed;

    
    
    public List<Source> smokeSources;
    public List<Source> targetSources;
    public float[] smoke0;
    public float[] smoke1;
    public float[] smoke0Blur;
    public float[] target;
    public float[] targetBlur;
    public boolean ifSmokeEpmty;
    
    public float DIFFUSSION = 1e-6f;
    public float VISICOSITY = 1e-6f;
  
    /**
     * initialize memory
     */
    public void setup() {
        elapsed = 0;
        N = Nval.getValue();        
        dx = 1.0f / N; // we choose the domain size here to be 1 unit square!
        
        int np2s = (N+2)*(N+2);
        U0 = new float[2][np2s];
        U1 = new float[2][np2s];
   
        smoke0 = new float[np2s];
        smoke0Blur = new float[np2s];
        smoke1 = new float[np2s];
        target = new float[np2s];
        targetBlur = new float[np2s];
        ifSmokeEpmty = true;  
        smokeSources = new LinkedList<Source>();
        targetSources = new LinkedList<Source>();
    }

    public int IX( int i, int j ) {
        return i*(N+2) + j;
    }
    
   
    public void setBoundary( int b, float[] x ) {
        int i;
        for ( i=1 ; i<=N; i++ ) {
            x[IX(0 ,i)]  = b==1  ? -x[IX(1,i)] : x[IX(1,i)];
            x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
            x[IX(i,0 )]  = b==2  ? -x[IX(i,1)] : x[IX(i,1)];
            x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];            
        }
        x[IX(0 ,0 )] = 0.5f*(x[IX(1,0 )]+x[IX(0 ,1)]);
        x[IX(0 ,N+1)] = 0.5f*(x[IX(1,N+1)]+x[IX(0 ,N )]);
        x[IX(N+1,0 )] = 0.5f*(x[IX(N,0 )]+x[IX(N+1,1)]);
        x[IX(N+1,N+1)] = 0.5f*(x[IX(N,N+1)]+x[IX(N+1,N )]);
    }

    
    
    public void getVelocity( Tuple2f x, Tuple2f vel ) {
        getVelocity( x, U0, vel );
    }
 
    private void getVelocity( Tuple2f x, float[][] U, Tuple2f vel ) {
        vel.x = interpolate( x, U[0] );
        vel.y = interpolate( x, U[1] );
    }
    
   
    public float interpolate( Tuple2f position, float[] s ) {
    	
    	// TODO: Objective 1: implement bilinear interpolation (try to make this code fast!)
 
    	int i0, j0, i1, j1;
    	float x_x1, x2_x, y_y1, y2_y;
    	
    	float x = position.x;
    	float y = position.y;
    	
    	if (x<0.5*dx) x=(float)0.5*dx; if (x>=(N+1.5)*dx) x=(float)((N+1.48)*dx); 
    	if (y<0.5*dx) y=(float)0.5*dx; if (y>=(N+1.5)*dx) y=(float)((N+1.48)*dx); 
    	    
    	i0=(int)(x/dx); 
    	if(x-i0*dx < 0.5*dx) {
    		i0 = i0-1;
    	}
    	i1=i0+1;
     	
    	j0=(int)(y/dx); 
    	if(y-j0*dx < 0.5*dx) {
    		j0 = j0-1;
    	}
    	j1=j0+1;
    	
    	x_x1 = (float) (x-((i0+1)*dx-0.5*dx)); 
    	x2_x = dx-x_x1; 
    	y_y1 = (float) (y-((j0+1)*dx-0.5*dx)); 
    	y2_y = dx-y_y1; 

    	float result = 0;
    	result += s[IX(i0,j0)] * x2_x *y2_y /(dx*dx);
    	result += s[IX(i1,j0)] * x_x1 *y2_y /(dx*dx);
    	result += s[IX(i0,j1)] * x2_x *y_y1 /(dx*dx);
    	result += s[IX(i1,j1)] * x_x1 *y_y1 /(dx*dx);
    	
    	return result;
    }
        
   
    public void traceParticle( Point2f x0, float h, Point2f x1 ) {
        traceParticle( x0, U0, h, x1 );        
    }
    
  
    private void traceParticle( Point2f x0, float[][] U, float h, Point2f x1 ) {
    
        Point2f u_x0 = new Point2f();
    	getVelocity(x0,U,u_x0);
    	x1.x  = x0.x + h * u_x0.x;
    	x1.y  = x0.y + h * u_x0.y;
    }
        

    public void transport( float[] s1, float[] s0, float[][] U, float dt ) {
        
    	Point2f prev = new Point2f();
    	Point2f current = new Point2f();
    	
    	
    	for ( int i=1 ; i<=N ; i++ ) {
    		for ( int j=1 ; j<=N ; j++ ) {
    			current.x = (float) ((i+0.5)*dx);
    			current.y = (float) ((j+0.5)*dx);
    			traceParticle(current, U, -dt, prev);   			
    			 	    
    			s1[IX(i,j)] = interpolate( prev, s0 );
    		}
    	}
    	
    }
    
    /**
     * Does the Poisson solve to make sure that velocities U respect incompressible flow
     * @param U
     */
    float[] div = new float[(N+2)*(N+2)];
	float[] p = new float[(N+2)*(N+2)];
    private void project( float[][] U ) {    

    	int i, j, k;
    	float h = dx;
    	for ( i=0; i<(N+2)*(N+2) ; i++ ) {
    		div[i] = 0;
    		p[i] = 0;
    	}
    	
    	for ( i=1 ; i<=N ; i++ ) {
    		for ( j=1 ; j<=N ; j++ ) {
    			div[IX(i,j)] = (float)( -0.5*h*(U[0][IX(i+1,j)]-U[0][IX(i-1,j)]+
            	    	                        U[1][IX(i,j+1)]-U[1][IX(i,j-1)]) );
            	p[IX(i,j)] = 0;
    		}
    	}
    	setBoundary (0, div);
    	setBoundary (0, p);
    	
    	for ( k=0 ; k< iterations.getValue() ; k++ ) {
    		for ( i=1 ; i<=N ; i++ ) {
    	    	for ( j=1 ; j<=N ; j++ ) {
    	    		p[IX(i,j)] = (div[IX(i,j)]+p[IX(i-1,j)]+p[IX(i+1,j)]+
    	       	    	 p[IX(i,j-1)]+p[IX(i,j+1)])/4;
    	    	}
    	    	
    	    }    	
    	    setBoundary (0, p );
    	}
    	for ( i=1 ; i<=N ; i++ ) {
    	    for ( j=1 ; j<=N ; j++ ) {
    	    	U[0][IX(i,j)] -= 0.5*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
    	    	U[1][IX(i,j)] -= 0.5*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
    	    }
    	}
    	setBoundary (1, U[0]);
    	setBoundary (2, U[1]);
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////addSource for target field and smoke field///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    
    public void addSource_target( float[] S, float dt, Tuple2f position) {
    	addSource( S, dt, position, 40f );
    }
    public void addSource_smoke( float[] S, float dt, Tuple2f position) {
    	float a = smokeAmount.getFloatValue();
    	addSource( S, dt, position, a);
    }
 
    public void addSource( float[] S, float dt, Tuple2f position, float amount ) {

    	
    	
    	//Distributed amount
    	float a11, a12, a21, a22;
    	// weight
    	float w11, w12, w21, w22;
    	
    	int i0, j0, i1, j1;
    	float x_x1, x2_x, y_y1, y2_y;
    	
    	float x = position.x;
    	float y = position.y;
    	
    	
    	if (x<0.5*dx) x=(float)0.5*dx; if (x>=(N+1.5)*dx) x=(float)((N+1.48)*dx); 
    	if (y<0.5*dx) y=(float)0.5*dx; if (y>=(N+1.5)*dx) y=(float)((N+1.48)*dx); 
    	    
    	i0=(int)(x/dx); 
    	if(x-i0*dx < 0.5*dx) {
    		i0 = i0-1;
    	}
    	i1=i0+1;
     	
    	j0=(int)(y/dx); 
    	if(y-j0*dx < 0.5*dx) {
    		j0 = j0-1;
    	}
    	j1=j0+1;
    	
    	x_x1 = (float) (x-((i0+1)*dx-0.5*dx)); 
    	x2_x = dx-x_x1; 
    	y_y1 = (float) (y-((j0+1)*dx-0.5*dx)); 
    	y2_y = dx-y_y1;
    	
    	w11 = x2_x *y2_y /(dx*dx);
    	w21 = x_x1 *y2_y /(dx*dx);
    	w12 = x2_x *y_y1 /(dx*dx);
    	w22 = x_x1 *y_y1 /(dx*dx);
    	
    	// amount = a11*w11 + a21*w21 + a12*w12 + a22*w22
    	// a11/w11 = a21/w21  -->  a21 = w21*a11/w11
        float a = Math.abs(amount);
    	
        if(w11!=0) {
    		a11 = a/(w11 + w21*w21/w11 + w12*w12/w11 + w22*w22/w11);
            a21 = w21*a11/w11;
            a12 = w12*a11/w11;
            a22 = w22*a11/w11;
    	}else if(w21!=0) {
    		a21 = a/(w21 + w11*w11/w21 + w12*w12/w21 + w22*w22/w21);
            a11 = w11*a21/w21;
            a12 = w12*a21/w21;
            a22 = w22*a21/w21;
    	}else if(w12!=0) {
    		a12 = a/(w12 + w11*w11/w12 + w21*w21/w12 + w22*w22/w12);
            a11 = w11*a12/w12;
            a21 = w21*a12/w12;
            a22 = w22*a12/w12;
    	}else {
    		a22 = a/(w22 + w11*w11/w22 + w21*w21/w22 + w12*w12/w22);
            a11 = w11*a22/w22;
            a12 = w12*a22/w22;
            a21 = w21*a22/w22;
    	}
        
        if(amount<0) {
        	a11 = a11*-1;
        	a21 = a21*-1;
        	a12 = a12*-1;
        	a22 = a22*-1;
        }
        
        S[IX(i0,j0)] += a11*dt;
    	S[IX(i1,j0)] += a21*dt;
    	S[IX(i0,j1)] += a12*dt;
    	S[IX(i1,j1)] += a22*dt;
    	
    	
        
    }
    
   
 
    
    /** Worker variables for mouse interaction */
    private Point2f XVX = new Point2f();
    private Point2f Xprev = new Point2f();
    
   
    
    /**
     * Sets the mouse location in the fluid for doing dragging interaction
     * @param x0
     * @param x1
     */
    public void setMouseMotionPos( Point2f x0, Point2f x1 ) {
        Xprev.set( x0 );
        XVX.set( x1 );
    }
    
    ///////////////////////////////////////////////////////////////////////////
    ////////////Vorticity added for smoke simulation/////////////////////////// 
    ///////////////////////////////////////////////////////////////////////////
    
      
    public void confineVorticity(float[] U_x, float[] U_y,float dt, float VORTICITY) {
    	
    	 // compute |w|, the curl, at each position in the velocity field
        float[] w = new float[(N+2)*(N+2)];
        
        for (int x=1; x<=N; x++) {
            for (int y=1; y<=N; y++) {
                w[IX(x,y)] = (float) ((U_y[IX(x+1,y)] - U_y[IX(x-1,y)] - 
                		U_x[IX(x,y+1)] + U_x[IX(x,y-1)]) / 2.0);              
            }
        }
        setBoundary(0,w);
        

        float dw_dy, dw_dx, norm;
        float[] fx_conf = new float[(N+2)*(N+2)];
        float[] fy_conf = new float[(N+2)*(N+2)];
        
        for (int y=1; y<=N; y++) {
            for (int x=1; x<=N; x++) {
            	dw_dx = (float) ((Math.abs(w[IX(x+1,y)]) - Math.abs(w[IX(x-1,y)])) / 2.0);
                dw_dy = (float) ((Math.abs(w[IX(x,y+1)]) - Math.abs(w[IX(x,y-1)])) / 2.0);
                
                norm = (float) Math.sqrt(dw_dy * dw_dy + dw_dx * dw_dx);
                dw_dx /= norm + 1e-5;
                dw_dy /= norm + 1e-5;
                
                fx_conf[IX(x,y)] = VORTICITY * dw_dy * -w[IX(x,y)];
                fy_conf[IX(x,y)] = VORTICITY * dw_dx * w[IX(x,y)];
               
            }
        }
        setBoundary(0,fy_conf);
        setBoundary(0,fx_conf);

        for (int x=1; x<=N; x++) {
            for (int y=1; y<=N; y++) {
            	U_x[IX(x,y)] += fx_conf[IX(x,y)] * dt;
                U_y[IX(x,y)] += fy_conf[IX(x,y)] * dt;             
            }
        }
        setBoundary(1,U_x);
        setBoundary(2,U_y);
    }
 
    
    ///////////////////////////////////////////////////////////////
    ///////////////////////Target driven///////////////////////////
    ///////////////////////////////////////////////////////////////
    
    //
    // Compute driving force
    // Vf is drving force parameter from the paper.
    //
    public void drive_force(float[] U1_x, float[] U1_y, float[] U0_x, float[] U0_y, float[] p, float[] target_p, float Vf,float dt) {
        
    	float Ap, Ap_star, Dp_star;
        for (int y = 1; y <=N; ++y) {
            for (int x = 1; x <=N; ++x) {
                // Vy
                Ap = (p[IX(x,y)] + p[IX(x,y+1)]) / 2.0f;
                Ap_star = (target_p[IX(x,y)] + target_p[IX(x,y+1)]) / 2.0f;
                Dp_star = target_p[IX(x,y+1)] - target_p[IX(x,y)];
                if (Ap_star == 0.0f) {
                    U1_y[IX(x,y)] = U0_y[IX(x,y)];
                } else {
                    U1_y[IX(x,y)] = U0_y[IX(x,y)] + dt * Vf * (Ap * Dp_star / Ap_star);
                }

                // Vx
                Ap = (p[IX(x,y)] + p[IX(x+1,y)]) / 2.0f;
                Ap_star = (target_p[IX(x,y)] + target_p[IX(x+1,y)]) / 2.0f;
                Dp_star = target_p[IX(x+1,y)] - target_p[IX(x,y)];
                if (Ap_star == 0.0f) {
                    U1_x[IX(x,y)] = U0_x[IX(x,y)];
                } else {
                    U1_x[IX(x,y)] = U0_x[IX(x,y)] + dt * Vf * (Ap * Dp_star / Ap_star);
                }
            }
        }
        setBoundary(1,U1_x);
        setBoundary(2,U1_y);
    }
    //
    // Momentum attenuation where vd is rate of momentum attenuation.
    //
    public void attenuate(float[] x, float[] x0, float vd, float dt) {
        // initial solution
    	for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {               
                x[IX(i,j)] =  x0[IX(i,j)];
            }
        }
        
        // solve by Gauss iterations
        int currIteration = 0;
        while (currIteration < iterations.getValue()) {
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {                 
                    x[IX(i,j)] /= (dt * vd + 1);
                }
            }
            currIteration ++;
        }
    }

    //
    // Smoke gathering 
    // Vg is the rate at which ¦Ñ is gathered toward ¦Ñ* from the paper.
    //
   
    public void gather( float[] S, float[] S0, float[] tg, float[] blr_tg, float vg, float dt) {
        
        // initial solution
    	for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {             
                S[IX(i,j)] =  S0[IX(i,j)];
            }
        }
        // solve by Gauss iterations
        int curr = 0;
        while (curr < iterations.getValue()) {
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    float rhs_x = 0;
                    float lhs_x = 0;
                  
                   // compute scales 
                    rhs_x += ( S[IX(i-1,j)]-tg[IX(i-1,j)]+tg[IX(i,j)] ) * blr_tg[IX(i,j)];
                    lhs_x += blr_tg[IX(i,j)];
             
                    rhs_x += (  S[IX(i+1,j)]-tg[IX(i+1,j)]+tg[IX(i,j)] ) * blr_tg[IX(i+1,j)] * S[IX(i+1,j)];
                    lhs_x += blr_tg[IX(i+1,j)] * S[IX(i+1,j)];
                
                    rhs_x += ( S[IX(i,j-1)]-tg[IX(i,j-1)]+tg[IX(i,j)] ) * blr_tg[IX(i,j)];
                    lhs_x += blr_tg[IX(i,j)];
                
                    rhs_x += ( S[IX(i,j+1)]-tg[IX(i,j+1)]+tg[IX(i,j)] )* blr_tg[IX(i,j+1)] * S[IX(i,j+1)];
                    lhs_x += blr_tg[IX(i,j+1)] * S[IX(i,j+1)];   
                    
                    // update solution
                    S[IX(i,j)] =  (float) (S0[IX(i,j)]*dx*dx + vg * dt * rhs_x);
                    S[IX(i,j)] /= ( dx*dx + vg*dt*lhs_x );
                }
            }
            curr ++;
        }
    }

    // A simple but stable Gaussian blur for Smoke field
    public void gaussian_blur1(float[] S, float[] S0) {
        for (int y = 0; y < N+2; y++) {
            for (int x = 0; x < N+2; x++) {
                if (y > 0 && x > 0)
                    S[IX(x, y)] += 1.0f * S0[IX(x-1,y-1)];
                if (y > 0)
                    S[IX(x, y)] += 2.0f * S0[IX(x,y-1)];
                if (y > 0 && x < N+1)
                    S[IX(x, y)] += 1.0f * S0[IX(x+1,y-1)];
                if (x > 0)
                    S[IX(x, y)] += 2.0f * S0[IX(x-1,y)];
                S[IX(x, y)]     += 4.0f * S0[IX(x,y)];
                if (x < N+1)
                    S[IX(x, y)] += 2.0f * S0[IX(x+1,y)];
                if (y < N+1 && x > 0)
                    S[IX(x, y)] += 1.0f * S0[IX(x-1,y+1)];
                if (y < N+1)
                    S[IX(x, y)] += 2.0f * S0[IX(x,y+1)];
                if (y < N+1 && x < N+1)
                    S[IX(x, y)] += 1.0f * S0[IX(x+1,y+1)];
                S[IX(x, y)] /= 16.0f;
            }
        }
    }
   
     // A complex version of Gaussian blur for target field
    public void gaussian_blur2(float[] d, float[] d0, float m_sigma) {
       
        float[][] d_arr = new float[N+2][N+2];
                
        for (int i = 0; i < (N+2); i++) {   
        	for (int j = 0; j < (N+2); j++) {
        		d_arr[i][j] = d0[IX(i,j)];  
        	}                      
        }
             
        GaussFilter2d(d_arr,  N+2,  N+2, m_sigma, 3);
             
        for (int i = 0; i < (N+2); i++) {   
        	for (int j = 0; j < (N+2); j++) {
        		d[IX(i,j)]=d_arr[i][j] ;  
        	}                      
        }
    }
    void GaussFilter2d(float[][] image, int width, int height,
            float sigma, int numsteps) 
    {
        
        double lambda, dnu;
        float nu, boundaryscale, postscale;       
        int i,j, x, y;
        int step;
               
        lambda = (sigma*sigma)/(2.0*numsteps);
        dnu = (1.0 + 2.0*lambda - Math.sqrt(1.0 + 4.0*lambda))/(2.0*lambda);
        nu = (float)dnu;
        boundaryscale = (float)(1.0/(1.0 - dnu));
        postscale = (float)(Math.pow(dnu/lambda,2*numsteps));
        
        /* Filter horizontally along each row */
        for(y = 0; y < height; y++) {
            for(step = 0; step < numsteps; step++) {
            	               
                image[y][0] *= boundaryscale;
                
                /* Filter rightwards */
                for(x = 1; x < width; x++)
                    image[y][x] += nu*image[y][x - 1];
                
                image[y][width-1] *= boundaryscale;
                
                /* Filter leftwards */
                for(x = width-1; x > 0; x--)
                	image[y][x - 1] += nu*image[y][x];
            }
        }
        
        /* Filter vertically along each column */
        for(x = 0; x < width; x++) {
            for(step = 0; step < numsteps; step++) {
            	
                image[0][x] *= boundaryscale;
                
                /* Filter downwards */
                for(i = 1; i < height; i ++)
                	image[i][x] += nu*image[i - 1][x];
                
                image[height - 1][x] *= boundaryscale;
                
                /* Filter upwards */
                for(i = height-1; i > 0; i --)
                    image[i - 1][x] += nu*image[i][x];
            }
        }
        
        for(i = 0; i < height; i++) {
        	for(j = 0; j < width; j++) {
        		 image[i][j] *= postscale;       	        
        	}
        }
    
    }

    
    ///////////////////////////////////////////////////////////
    ////////////////// Step ///////////////////////////////////
    //////////////////////////////////////////////////////////
    public void velocityStep_target( float dt ) {
    	
    	float[][] tmp;
    	
    	// add driving force
    	gaussian_blur1(smoke0Blur, smoke0);   	
    	float Vf = driving_force.getFloatValue();
    	drive_force(U1[0], U1[1], U0[0], U0[1], smoke0Blur,targetBlur, Vf, dt);
    	tmp = U1; U1 = U0; U0 = tmp;
    	setBoundary( 1, U0[0] ); 
        setBoundary( 2, U0[1] );
        
        // add momentum attenuation
    	float Vd = attenuation.getFloatValue();
    	attenuate(U1[0],U0[0], Vd, dt);
        attenuate(U1[1],U0[1],Vd, dt);   
        tmp = U1; U1 = U0; U0 = tmp;
        setBoundary( 1, U0[0] ); 
        setBoundary( 2, U0[1] );
      
        // self transport
        transport( U1[0], U0[0], U1, dt );
        transport( U1[1], U0[1], U1, dt );
        tmp = U1; U1 = U0; U0 = tmp;
        setBoundary( 1, U0[0] ); 
        setBoundary( 2, U0[1] );
        
        // confine vorticity
        float vort = vorticity.getFloatValue();
    	confineVorticity(U0[0], U0[1], dt, vort); 
    	
    	// projection
        project ( U0 );
        
    }
  
    
    public void scalarStep_target(float dt) {
    	float[] tmpt;  
    	
    	// Velocity advect	
	    transport( smoke1, smoke0, U0, dt );
	    tmpt = smoke1; smoke1 = smoke0; smoke0 = tmpt;
	    setBoundary(0, smoke0); 	
    	
	    // Gathering
	    if(!ifSmokeEpmty) {
	    	float Vg = gather_rate.getFloatValue();
	    	gather(smoke1, smoke0, target, targetBlur, Vg, dt);
	    	tmpt = smoke1; smoke1 = smoke0; smoke0 = tmpt;
		    setBoundary(0, smoke0);
	    }
       
    }
    
    public void step() {
        float dt = timeStepSize.getFloatValue();
        velocityStep_target(dt);
        scalarStep_target(dt);                
        elapsed += dt;
 
        
    }
    public void step_drawing() {
        gaussian_blur2(targetBlur, target,1);       
    }
        
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////Parameters can be controlled////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    private DoubleParameter vorticity = new DoubleParameter( "vorticity", 1, 1e-3, 50 );
    private DoubleParameter gather_rate = new DoubleParameter( "Gather Rate", 7e-8, 1e-8, 1e-5f );
    private DoubleParameter driving_force = new DoubleParameter( "Driving Force", 1.3e-2, 1e-3, 10f);
    private DoubleParameter attenuation = new DoubleParameter( "Momentum Attenuation", 3.8e-2, 1e-2, 1);
    private DoubleParameter smokeAmount = new DoubleParameter( "Add Amount", 1800, 1000, 15000);
    
    private IntParameter iterations = new IntParameter( "GS iterations", 30, 0, 1000 );    
    public DoubleParameter timeStepSize = new DoubleParameter( "time step size", 0.1, 0.001, 1 );
    
  
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();

        VerticalFlowPanel vfp1 = new VerticalFlowPanel();
        vfp1.setBorder(new TitledBorder("Fluid properties"));
       
        vfp1.add( vorticity.getSliderControls(true) );
        vfp1.add( gather_rate.getSliderControls(true) );       
        vfp1.add( driving_force.getSliderControls(true) );
        vfp1.add( attenuation.getSliderControls(true) );   
        vfp1.add( smokeAmount.getSliderControls(true) );      
        vfp.add( vfp1.getPanel() );

        VerticalFlowPanel vfp2 = new VerticalFlowPanel();
        vfp2.setBorder(new TitledBorder("Fluid solver properties"));    
        vfp2.add( timeStepSize.getSliderControls(true ) ); 
        vfp2.add( iterations.getSliderControls() );
        vfp.add( vfp2.getPanel() );
 
        		
        return vfp.getPanel();
    }
}