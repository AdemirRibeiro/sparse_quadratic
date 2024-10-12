function qt_outer()
  sv=true; nf=0; lw=2; szft=10;
  cor0=RGB(1,1,1); cor1=RGB(.1,.1,.9);
  cor2=RGB(.8,.1,.2); cor3=RGB(.25,.25,.25); cor4=RGB(1,0,1);
  s=4;
  qt=[];
  for n=[12]
    j=0:n;
    j1=ceil(Int,(n+1)/2);
    qt=factorial.(n,n.-j)./factorial.(j);
    aux=(scatter(j,qt,color=cor1,xtickfont=font(szft),ytickfont=font(szft),
         xtick=(0:n),ytick=([qt[1:j1];2^n]),legend=false);
         # annotate!(n,qt[j1],2^n,szft)
         )
    display(aux)
    if sv nf=nf+1; savefig("/home/ademir/Ademir/artigos/Mael/fig/qt_outer$nf.pdf"); end

    # aux=(scatter!([s],[qt[s+1]],color=cor2);)
    # display(aux)
    # if sv nf=nf+1; savefig("/home/ademir/Ademir/artigos/Mael/fig/qt_outer$nf.pdf"); end

  end
end

qt_outer()
