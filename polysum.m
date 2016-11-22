function s=polysum(p,q)
maxlen=max(lenght(p)), linght(q));
p1=zeros(1, maxlen);
q1=zeros(1, maxlen);
p1(maxlen-length(p)+1:maxlen)=p;
q1(maxlen-length(q)+1:maxlen)=q;
s=p1+q1;