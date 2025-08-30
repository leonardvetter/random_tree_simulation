# 1. Setup

1. Clone repo
2. conda env create --file=environment.yml
3. open demo_sim_random_tree.ipynb

# 2. Theory of simulating random plane trees

This serves as an introduction to simulating properties of size-conditioned BGW trees, via the simulation of certain random walks. There is a one-to-one correspondence between plane trees and integer valued walks. The parameter for the simulation is the offspring distribution $\mu$ and the size (i.e. the number of vertices) of a tree.
        
Thus, we can simulate random trees by simulating certain random walks and translating them to trees afterwards. With this in mind let $\mu$ be an offspring distribution and $X$ be a random variable with law given by $\mathbb P(X=k)=\mu_{k+1}$ for $k\in \{-1,0,1,2,\dots\}$. Let $(X_i)$ be i.i.d. distributed as $X$, define $W_n$ to be the sum of $X_i$ from 1 to $n$ and $\zeta= \inf ( k\geq 0:W_k=-1 )$. Recall that the Lukasiewicz path of a $\textrm{BGW}_\mu$ tree has the same distribution as the random walk, stopped at $\zeta$.
        
We study $T_n$, which is a $\textrm{BGW}_\mu$ tree under the conditional probability that the tree has $n$ vertices. This conditioning translates to studying $(W_0,\dots,W_n)$ under the conditional probability that the first time, the random walk hits $-1$ is at time $n$. This means that the random walk is conditioned to stay non-negative for $1\leq i<n$ and $W_n=-1$. It is difficult to simulate under the excursion condition since the probability of this event is small.
        
Fortunately, we can weaken the conditioning by using the Vervaat-transform. The Vervaat-transform is basically a cyclic shift of the increments of the random walk by the index of the first time the random walk hits its global minima. More formally, let $n\in \mathbb{N}$, $(x_1,\dots,x_n)\in \mathbb{Z}^n$ and define the corresponding (deterministic) walk

       
$$w_i = \sum_{j=1}^i x_j, \hspace{1cm} 1\leq i\leq j.$$

Furthermore, define the first time at which the random walk reaches its overall minimum by

$$M_n := \min \{0\leq i\leq n: w_i = \min \{ w_j:0\leq j \leq n\} \}.$$

The Vervaat transform is a cyclic shift by $M_n$ of the increments of the random walk. The cyclic shift of $x = (x_1,\dots,x_n)\in \mathbb Z^n$ with $M_n$ is defined to be 
         $$x^{(M_n)}:=(x_{1+M_n} ,x_{2+M_n},\dots,x_{n+M_n})$$
        with the addition in the index modulo $n$. The Vervaat transform $\mathcal{V}(w)$ is the corresponding random walk of the shifted increments, i.e. $\mathcal{V}(w) = \sum_{i=1}^n x^{(M_n)}_i$. The main reason why the Vervaat transform is useful for us is the following well known fact, see e.g. \cites{Pit06}.
        
The law of $(X_1,\dots,X_n)$ under $\mathbb P(\cdot\mid \zeta=n)$ is equal to the law of $\mathcal{V}\big((X_1,\dots,X_n)\big)$ under $\mathbb P(\cdot\mid W_n=-1)$.
        
Hence, instead of simulating under the conditional probability $\mathbb P(\cdot\mid \zeta=n)$, we can simulate under the conditional probability $\mathbb P(\cdot\mid W_n=-1)$ and apply the Vervaat transform to the result. We define the bridge $B^{(n)}:= (B_i^{(n)}:0\leq i\leq n)$ to be distributed as $(W_i:0\leq i \leq n)$ under the conditional probability $\mathbb P(\cdot \mid W_n = -1)$. For $i \in \{0,\dots,n-1\}$ denote by $b_i^{(n)}:=B_{i+1}^{(n)}-B_i^{(n)}$ the increments of the bridge. The algorithm for sampling the random excursion is basically the following rejection sampling.

1. Sample the increments $(X_1,\dots,X_n)$ with law given by $\mathbb P(X=k)=\mu_{k+1}$ for $k\in \{-1,0,1,2,\dots\}$ independent of each other.

2. Reject if $\sum_{k=1}^n X_i \neq -1$ and return to 1.

3. Apply the Vervaat-transform $\mathcal{V}((X_1,\dots,X_n))$.

For the first step, [Devroye12](https://www.researchgate.net/publication/220617715_Simulating_Size-constrained_GaltonWatson_Trees) proposed an efficient sampling method. The idea is to sample for each $i\in \{-1,0,1,2,\dots\}$ the number of jumps of size $i$, instead of sampling the jumps directly. Then, one creates a vector of the jumps in some order and applies a random permutation.
        
For this method, one samples a multinomial random vector $(N_{-1},N_0,N_1,N_2,\dots)$ with parameters $(n,\mu_0,\mu_1,\mu_2,\mu_3,\dots)$. For $i \geq -1$ the random variable $N_i$ is the number of jumps of size $i$. We can calculate the sum of the increments by 


$$W_n = \sum_{i=-1}^\infty i N_i.$$


There is one small issue remaining, since the offspring distribution $(\mu_k)_{k\geq 0}$ has to be truncated for the simulation. We can prove easily that if we simulate $\mathcal T_n$ under the truncated offspring distribution, the simulation is still is exact. More precisely, let $k\geq n$ and define the \textbf{truncated offspring distribution} $\mu^{(k)} = (\mu^{(k)}_i)_{i\geq 0}$ by


$$\mu^{(k)}_i:= \frac{\mu_i}{\mu([0,k])}, \textrm{ if } i\in \{0,\dots,k\} \textrm{ and } \mu^{(k)}_i:= 0 \textrm{ if } i>k.$$

In the demo notebook, we use an example offspring distribution that is subcritical (expected number of children of one individual = m < 1) and belongs to the domain of attraction of a stable law of index $\alpha = 1$. When we condition on the size of the tree being equal to $n$, an interesting phenomenon arises, as $n \rightarrow \infty$: One vertex with macroscopic degree $(1-m)n$ emerges. (For more details see: (https://arxiv.org/abs/2503.07530))



