Lam=NaN*ones(I,J,180);

half_dt=2;
dt=2*half_dt+1;
win_delta=2;
w=zeros(dt,(2*win_delta+1)^2);
for start_pos=3:139
for i=win_delta+1:I-win_delta
    for j=win_delta+1:J-win_delta
        for k=1:dt
            s=squeeze(S(i-win_delta:i+win_delta,j-win_delta:j+win_delta,start_pos-1+k));
            w(k,:)=s(:)';
        end
        mw=mean(w,2);
        ww=w-repmat(mw,1,(2*win_delta+1)^2);
        A=ww*ww';
        e1=A*A*ones(dt,1);
        e1=e1/sqrt(e1'*e1);
        r=ww-e1*e1'*ww;
        Lam(i,j,start_pos+half_dt)=mean(sqrt(sum(r.*r,2)));
    end
end
end



