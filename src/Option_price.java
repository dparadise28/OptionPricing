/*
 * Simple American options pricing model that assumes GBM 
 * -------------------------------------------------------------
 * pricing is done on binomial trees by recursively looking back
 * and finding expected payoffs until last node is reached. The 
 * result is the expected payoff of the final node
 */

public class Option_price {
	public static void main(String[] args)
	{
		long startTime = System.nanoTime();
		double[] Put_Call_price=americanPC(1,100,100,.05,.2,.0,20000);
		
		System.out.print("put price: "+Put_Call_price[0]+"\n"+
						 "call price: "+Put_Call_price[1]+"\n"+
						 "Time elapsed (nano secs): "+(System.nanoTime()-startTime)*1.0/1000000000);
		
	}
	
	//computes the price of a put and a call using a binomial tree of size n
	private static double[] americanPC(double t, double S, double K, double r, double sigma, double q, int n) 
	{
	    double dt=t*1.0/n, u=Math.exp(sigma*Math.sqrt(t*1.0/n)), d=Math.exp(-sigma*Math.sqrt(t*1.0/n));
	    	   /*
	    	    * n=number if time periods
	    	    *                             
	    	    *                               |    S0*u^n 
	    	    *                         ...   |
	    	    *                               |    S0*u^(n-1)*d
	    	    *                         ...   |
	    	    *                 S0*u^2        |        .
	    	    *        u*S0             ...   |        .
	    	    * s0              S0*u*d        |        .
	    	    *        d*S0             ...   |
	    	    *                 S0*d^2        |
	    	    *                         ...   |    S0*u*d^(n-1)
	    	    *                               |    
	    	    *                         ...   |    S0*d^n
	    	    *                               |
	    	    *                               |
	    	    * t=0     t=1      t=2    ...        t=n      
	    	    * 
	    	    * dt=time elapsed between tree layers
	    	    * u=multiplicative up factor
	    	    * d=multiplicative down factor
	    	    */
	    
	    double p=(Math.exp((r-q)*dt)-d)/(u-d);  //probability of the underlying price increase
	    double[] p_tn=new double[n+1];          //array used to store possible prices at time n
	    
	    //compute all possible prices at time n (maturity)
	    for(int i=n; i>=0; i--)
	    	p_tn[i]=Math.max(K-S*Math.pow(u,n-i)*Math.pow(d,i), 0);
	    
	    //recursively loop through all paths of the tree and price the option at every node by finding 
	    //expected payoff and early exercise
	    for(int j=n-1; j>=0; j--)
	        for(int i=0; i<=j; i++) 
	        	p_tn[i]=Math.max(Math.exp(-(r-q)*dt)*(p*p_tn[i]+(1-p)*p_tn[i+1]), K-S*Math.pow(u,n-i)*Math.pow(d,i)); 
	    		// p*p_tn[i]+(1-p)*p_tn[i+1]=price of the put
	            // K-S*Math.pow(u,n-i)*Math.pow(d,i)=early exercise value of put
	    		// p_tn[i]=max(expected payoff, K-S_t); if no arbitrage is possible
	    
	    
	    //return put and call price (call price compute using put call parity)
	    return P_C_Parity(p_tn[0],S,K,r,t,q); 
	}
	
	private static double[] P_C_Parity(double put, double S, double k, double r, double t, double q)
	{
		//put call parity (with dividends)
		double put_call[]={put, put-k*Math.exp(-r*t)+S*Math.exp(-q*t)};
		
		return put_call;
	}
}
