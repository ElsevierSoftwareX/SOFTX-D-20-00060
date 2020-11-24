!bc

BV(1)=1/(2*b*c)*(8*b*c+(b-c-d)*abs04-(b+c-d)*abs02-b*abs03+c*abs03-d*abs03+b*abs01+c*abs01+d*abs01)

BV(2)=0

BV(3)=1/(6*b*c)*(8*b*c+(b-c-d)*abs04-(b+c-d)*abs02-b*abs03+c*abs03-d*abs03+b*abs01+c*abs01+d*abs01)

BV(4)=1/(6*b2*c)*((-2*b2+b*(c+d)+(c+d)**2)*abs04+(2*b2+b*(c-d)-(c-d)**2)*abs02-2*b2*abs03+b*c*abs03+&
c2*abs03-b*d*abs03-2*c*d*abs03+d2*abs03+2*b2*abs01+b*c*abs01-c2*abs01+b*d*abs01-2*c*d*abs01-d2*abs01)

BV(5)=0

BV(6)=1/(18*b2*c)*((-2*b2+b*(c+d)+(c+d)**2)*abs04+(2*b2+b*(c-d)-(c-d)**2)*abs02-2*b2*abs03+b*c*abs03+c2*abs03-b*d*abs03-2*c*d*abs03+d2*abs03+2*b2*abs01+b*c*abs01-&
c2*abs01+b*d*abs01-2*c*d*abs01-d2*abs01)

BV(7)=1/(12*b3*c)*(16*b3*c+(3*b3-b2*(c+d)-b*(c+d)**2-(c+d)**3)*abs04+(-3*b3+b*(c-d)**2-(c-d)**3+&
b2*(-c+d))*abs02-3*b3*abs03+b2*c*abs03+b*c2*abs03+c3*abs03-b2*d*abs03-2*b*c*d*abs03-3*c2*d*abs03+b*d2*abs03+3*c*d2*abs03-d3*abs03+3*b3*abs01+&
b2*c*abs01-b*c2*abs01+c3*abs01+b2*d*abs01-2*b*c*d*abs01+3*c2*d*abs01-b*d2*abs01+3*c*d2*abs01+d3*abs01)

BV(8)=0

BV(9)=1/(36*b3*c)*(16*b3*c+(3*b3-b2*(c+d)-b*(c+d)**2-(c+d)**3)*abs04+(-3*b3+b*(c-d)**2-(c-d)**3+b2*(-c+d))*abs02-3*b3*abs03+b2*c*abs03+b*c2*abs03+c3*abs03-b2*d*abs03-&
2*b*c*d*abs03-3*c2*d*abs03+b*d2*abs03+3*c*d2*abs03-d3*abs03+3*b3*abs01+b2*c*abs01-b*c2*abs01+c3*abs01+b2*d*abs01-2*b*c*d*abs01+3*c2*d*abs01-&
b*d2*abs01+3*c*d2*abs01+d3*abs01)

BV(10)=1/(6*b*c2)*((b2-2*c2+b*(c-2*d)-c*d+d2)*abs04-(b2-2*c2+c*d+d2-b*(c+2*d))*abs02+b2*abs03+b*c*abs03-&
2*c2*abs03+2*b*d*abs03+c*d*abs03+d2*abs03-b2*abs01+b*c*abs01+2*c2*abs01-2*b*d*abs01+c*d*abs01-d2*abs01)

BV(11)=0

BV(12)=1/(18*b*c2)*((b2-2*c2+b*(c-2*d)-c*d+d2)*abs04-(b2-2*c2+c*d+d2-b*(c+2*d))*abs02+b2*abs03+b*c*abs03-2*c2*abs03+2*b*d*abs03+c*d*abs03+d2*abs03-b2*abs01+b*c*abs01+&
2*c2*abs01-2*b*d*abs01+c*d*abs01-d2*abs01)

BV(13)=1/(24*b2*c2)*((-3*b3+(3*c-d)*(c+d)**2+b2*(-3*c+5*d)+b*(3*c2+2*c*d-d2))*abs04+(3*b3+&
(c-d)**2*(3*c+d)-b2*(3*c+5*d)+b*(-3*c2+2*c*d+d2))*abs02+3*b3*abs03+3*b2*c*abs03-3*b*c2*abs03-3*c3*abs03+5*b2*d*abs03+2*b*c*d*abs03+5*c2*d*abs03+&
b*d2*abs03-c*d2*abs03-d3*abs03-3*b3*abs01+3*b2*c*abs01+3*b*c2*abs01-3*c3*abs01-5*b2*d*abs01+2*b*c*d*abs01-5*c2*d*abs01-b*d2*abs01-c*d2*abs01+&
d3*abs01)

BV(14)=0

BV(15)=1/(72*b2*c2)*((-3*b3+(3*c-d)*(c+d)**2+b2*(-3*c+5*d)+b*(3*c2+2*c*d-d2))*abs04+(3*b3+(c-d)**2*(3*c+d)-b2*(3*c+&
5*d)+b*(-3*c2+2*c*d+d2))*abs02+3*b3*abs03+3*b2*c*abs03-3*b*c2*abs03-3*c3*abs03+5*b2*d*abs03+2*b*c*d*abs03+5*c2*d*abs03+b*d2*abs03-c*d2*abs03-&
d3*abs03-3*b3*abs01+3*b2*c*abs01+3*b*c2*abs01-3*c3*abs01-5*b2*d*abs01+2*b*c*d*abs01-5*c2*d*abs01-b*d2*abs01-c*d2*abs01+d3*abs01)

BV(16)=1/(60*b3*c2)*((6*b4+b3*(6*c-9*d)-b*(4*c-d)*(c+d)**2-(4*c-d)*(c+d)**3+b2*(-4*c2-3*c*d+d2))*abs04-(6*b4+b*(c-d)**2*(4*c+d)-(c-d)**3*(4*c+d)-3*b3*(2*c+3*d)+&
b2*(-4*c2+3*c*d+d2))*abs02+6*b4*abs03+6*b3*c*abs03-4*b2*c2*abs03-4*b*c3*abs03-4*c4*abs03+9*b3*d*abs03+3*b2*c*d*abs03+7*b*c2*d*abs03+11*c3*d*abs03+&
b2*d2*abs03-2*b*c*d2*abs03-9*c2*d2*abs03-b*d3*abs03+c*d3*abs03+d4*abs03-6*b4*abs01+6*b3*c*abs01+4*b2*c2*abs01-4*b*c3*abs01+4*c4*abs01-&
9*b3*d*abs01+3*b2*c*d*abs01-7*b*c2*d*abs01+11*c3*d*abs01-b2*d2*abs01-2*b*c*d2*abs01+9*c2*d2*abs01+b*d3*abs01+c*d3*abs01-d4*abs01)

BV(17)=0

BV(18)=-1/(12*b*c3)*(-16*b*c3-(b3-3*c3+b2*(c-3*d)-c2*d+c*d2-d3+b*(c2-2*c*d+3*d2))*abs04-(-b3-3*c3+c2*d+c*d2+d3+b2*(c+3*d)-b*(c2+2*c*d+3*d2))*abs02+&
b3*abs03+b2*c*abs03+b*c2*abs03-3*c3*abs03+3*b2*d*abs03+2*b*c*d*abs03+c2*d*abs03+3*b*d2*abs03+c*d2*abs03+d3*abs03-b3*abs01+b2*c*abs01-&
b*c2*abs01-3*c3*abs01-3*b2*d*abs01+2*b*c*d*abs01-c2*d*abs01-3*b*d2*abs01+c*d2*abs01-d3*abs01)

BV(19)=0

BV(20)=-1/(36*b*c3)*(-16*b*c3-(b3-3*c3+b2*(c-3*d)-c2*d+c*d2-d3+b*(c2-2*c*d+3*d2))*abs04-(-b3-3*c3+c2*d+c*d2+d3+b2*(c+3*d)-b*(c2+2*c*d+3*d2))*abs02+b3*abs03+b2*c*abs03+&
b*c2*abs03-3*c3*abs03+3*b2*d*abs03+2*b*c*d*abs03+c2*d*abs03+3*b*d2*abs03+c*d2*abs03+d3*abs03-b3*abs01+b2*c*abs01-b*c2*abs01-3*c3*abs01-&
3*b2*d*abs01+2*b*c*d*abs01-c2*d*abs01-3*b*d2*abs01+c*d2*abs01-d3*abs01)

BV(21)=1/(60*b2*c3)*((-4*b4+b3*(-4*c+11*d)+b2*(-4*c2+7*c*d-9*d2)+&
(c+d)**2*(6*c2-3*c*d+d2)+b*(6*c3+3*c2*d-2*c*d2+d3))*abs04+(4*b4-b3*(4*c+11*d)-(c-d)**2*(6*c2+3*c*d+d2)+b2*(4*c2+7*c*d+9*d2)+b*(6*c3-3*c2*d-&
2*c*d2-d3))*abs02-4*b4*abs03-4*b3*c*abs03-4*b2*c2*abs03+6*b*c3*abs03+6*c4*abs03-11*b3*d*abs03-7*b2*c*d*abs03-3*b*c2*d*abs03-9*c3*d*abs03-&
9*b2*d2*abs03-2*b*c*d2*abs03+c2*d2*abs03-b*d3*abs03+c*d3*abs03+d4*abs03+4*b4*abs01-4*b3*c*abs01+4*b2*c2*abs01+6*b*c3*abs01-6*c4*abs01+&
11*b3*d*abs01-7*b2*c*d*abs01+3*b*c2*d*abs01-9*c3*d*abs01+9*b2*d2*abs01-2*b*c*d2*abs01-c2*d2*abs01+b*d3*abs01+c*d3*abs01-d4*abs01)

BV(22)=0

BV(23)=1/(180*b3*c3)*(80*b3*c3+(10*b5+2*b4*(5*c-13*d)-b*(c+d)**2*(10*c2-4*c*d+d2)-(c+d)**3*(10*c2-4*c*d+d2)+b3*(10*c2-16*c*d+19*d2)-b2*(10*c3+6*c2*d-3*c*d2+d3))*abs04+&
(-10*b5+2*b4*(5*c+13*d)+b*(c-d)**2*(10*c2+4*c*d+d2)-(c-d)**3*(10*c2+4*c*d+d2)-b3*(10*c2+16*c*d+19*d2)+b2*(-10*c3+6*c2*d+3*c*d2+d3))*abs02-&
10*b5*abs03-10*b4*c*abs03-10*b3*c2*abs03+10*b2*c3*abs03+10*b*c4*abs03+10*c5*abs03-26*b4*d*abs03-16*b3*c*d*abs03-6*b2*c2*d*abs03-16*b*c3*d*abs03-&
26*c4*d*abs03-19*b3*d2*abs03-3*b2*c*d2*abs03+3*b*c2*d2*abs03+19*c3*d2*abs03-b2*d3*abs03+2*b*c*d3*abs03-c2*d3*abs03+b*d4*abs03-c*d4*abs03-&
d5*abs03+10*b5*abs01-10*b4*c*abs01+10*b3*c2*abs01+10*b2*c3*abs01-10*b*c4*abs01+10*c5*abs01+26*b4*d*abs01-16*b3*c*d*abs01+6*b2*c2*d*abs01-&
16*b*c3*d*abs01+26*c4*d*abs01+19*b3*d2*abs01-3*b2*c*d2*abs01-3*b*c2*d2*abs01+19*c3*d2*abs01+b2*d3*abs01+2*b*c*d3*abs01+c2*d3*abs01-b*d4*abs01-&
c*d4*abs01+d5*abs01)

