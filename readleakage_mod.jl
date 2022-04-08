function read_leakage(Δm, Δl, ℓ, ℓ′)


  rleaks1 = OffsetArray(zeros(Float64,2*Δm+1,2*Δl+1,2*ℓ+1),-Δm:Δm,-Δl:Δl,-ℓ:ℓ)
  rleaks2 = OffsetArray(zeros(Float64,2*Δm+1,2*Δl+1,2*ℓ′+1),-Δm:Δm,-Δl:Δl,-ℓ′:ℓ′)
  horleaks1 = OffsetArray(zeros(Float64,2*Δm+1,2*Δl+1,2*ℓ+1),-Δm:Δm,-Δl:Δl,-ℓ:ℓ)
  horleaks2 = OffsetArray(zeros(Float64,2*Δm+1,2*Δl+1,2*ℓ′+1),-Δm:Δm,-Δl:Δl,-ℓ′:ℓ′)

  dl_mat=6
  dm_mat=15

  for ii in 1:2
    if(ii==1) lp=ℓ;end;
    if(ii==2) lp=ℓ′;end;
    ind0 = Int(lp*(lp+1)/2)

    for mp in -lp:lp
     ind = ind0 + mp + 1
     for j = -Δl:Δl
      l = lp + j

      for i in -Δm:Δm

       m = mp + i
       if (abs(m) > l) continue;end;


	   csphase = 1
	   if ((m>0)*(m%2==1))
	   	  csphase *= -1
       end

	   if ((mp>0)*(mp%2==1))
	   	  csphase *= -1
	   end

	   L_ii = ci_mat[(abs(m)-abs(mp))+dm_mat+1,l-lp+dl_mat+1,ind]*csphase
       L_rr = cr_mat[(abs(m)-abs(mp))+dm_mat+1,l-lp+dl_mat+1,ind]*csphase

       if (ii==1) 
	       rleaks1[i,j,mp] = (L_rr + sign(m*mp)*L_ii)/2
       else
	       rleaks2[i,j,mp] = (L_rr + sign(m*mp)*L_ii)/2
       end

       L_ii = hi_mat[abs(m)-abs(mp)+dm_mat+1,l-lp+dl_mat+1,ind]*csphase
       L_rr = hr_mat[abs(m)-abs(mp)+dm_mat+1,l-lp+dl_mat+1,ind]*csphase

       if (ii==1) 
	       horleaks1[i,j,mp] = (L_rr + sign(m*mp)*L_ii)/2
       else
	       horleaks2[i,j,mp] = (L_rr + sign(m*mp)*L_ii)/2
       end

      end
     end
    end
  end

  return rleaks1, horleaks1, rleaks2, horleaks2
end

#-------------------------------------------------------


