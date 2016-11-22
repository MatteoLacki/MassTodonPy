y <- c(100, 200, 50, 100, 250)

p = c(.3, .7, .7+.2, 0, 0 )
q = c(0, 0, .7+.2, .2, .8 )
P = matrix( c(p,q), ncol=2 )

GrammVector <- - t(P) %*% y
GrammMatrix <- t(P) %*% P

kernlab::ipop(	
	c 	= GrammVector, 
	H 	= GrammMatrix, 
	l 	= rep.int(0, ncol(P)), 
	u	= rep.int(100000, ncol(P)),
	A 	= t( matrix(rep.int(1, ncol(P))) ), 
	b 	= 0, 
	r 	= 10000000
)