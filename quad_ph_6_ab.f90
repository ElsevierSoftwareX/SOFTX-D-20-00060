! ab

BV(1)=1/(4*a*b)*(8*a*b+(a-b-c)*abs01-(a+b-c)*abs02-a*abs03+b*abs03-c*abs03+a*abs04+b*abs04+c*abs04)
BV(2)=1/(12*a2*b)*((-2*a2+a*(b+c)+(b+c)**2)*abs01+(2*a2+a*(b-c)-(b-c)**2)*abs02-2*a2*abs03+&
a*b*abs03+b2*abs03-a*c*abs03-2*b*c*abs03+c2*abs03+2*a2*abs04+a*b*abs04-b2*abs04+a*c*abs04-2*b*c*abs04-c2*abs04)
BV(3)=1/(24*a3*b)*(16*a3*b+(3*a3-a2*(b+c)-a*(b+c)**2-(b+c)**3)*abs01+(-3*a3+a*(b-c)**2-(b-c)**3+&
a2*(-b+c))*abs02-3*a3*abs03+a2*b*abs03+a*b2*abs03+b3*abs03-a2*c*abs03-2*a*b*c*abs03-3*b2*c*abs03+a*c2*abs03+3*b*c2*abs03-c3*abs03+3*a3*abs04+&
a2*b*abs04-a*b2*abs04+b3*abs04+a2*c*abs04-2*a*b*c*abs04+3*b2*c*abs04-a*c2*abs04+3*b*c2*abs04+c3*abs04)
BV(4)=1/(12*a*b2)*((a2-2*b2+a*(b-2*c)-b*c+c2)*abs01-(a2-2*b2+b*c+c2-a*(b+2*c))*abs02+a2*abs03+a*b*abs03-2*b2*abs03+2*a*c*abs03+b*c*abs03+&
c2*abs03-a2*abs04+a*b*abs04+2*b2*abs04-2*a*c*abs04+b*c*abs04-c2*abs04)
BV(5)=1/(48*a2*b2)*((-3*a3+(3*b-c)*(b+c)**2+a2*(-3*b+5*c)+a*(3*b2+2*b*c-c2))*abs01+(3*a3+(b-c)**2*(3*b+&
c)-a2*(3*b+5*c)+a*(-3*b2+2*b*c+c2))*abs02+3*a3*abs03+3*a2*b*abs03-3*a*b2*abs03-3*b3*abs03+5*a2*c*abs03+2*a*b*c*abs03+5*b2*c*abs03+a*c2*abs03-&
b*c2*abs03-c3*abs03-3*a3*abs04+3*a2*b*abs04+3*a*b2*abs04-3*b3*abs04-5*a2*c*abs04+2*a*b*c*abs04-5*b2*c*abs04-a*c2*abs04-b*c2*abs04+c3*abs04)
BV(6)=-1/(24*a*b3)*(-16*a*b3-(a3-3*b3+a2*(b-3*c)-b2*c+b*c2-c3+a*(b2-2*b*c+3*c2))*abs01-(-a3-3*b3+b2*c+b*c2+c3+a2*(b+3*c)-a*(b2+2*b*c+3*c2))*abs02+&
a3*abs03+a2*b*abs03+a*b2*abs03-3*b3*abs03+3*a2*c*abs03+2*a*b*c*abs03+b2*c*abs03+3*a*c2*abs03+b*c2*abs03+c3*abs03-a3*abs04+a2*b*abs04-&
a*b2*abs04-3*b3*abs04-3*a2*c*abs04+2*a*b*c*abs04-b2*c*abs04-3*a*c2*abs04+b*c2*abs04-c3*abs04)
