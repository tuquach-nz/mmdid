/*
Run this code to simulation
Where: mm is the ATT estimated by Multilevel Modeling approach
	   cs is the ATT estimated by Callaway and Sant'Anna approach

run "mmdid"

N = 100
simulate mm5=r(att5) mm6=r(att6) mm7=r(att7) cs5=r(cs5) cs6=r(cs6) cs7=r(cs7), rep(10000): mmdid, n(10)

N=500
simulate mm5=r(att5) mm6=r(att6) mm7=r(att7) cs5=r(cs5) cs6=r(cs6) cs7=r(cs7), rep(10000): mmdid, n(50)

N=1000
simulate mm5=r(att5) mm6=r(att6) mm7=r(att7) cs5=r(cs5) cs6=r(cs6) cs7=r(cs7), rep(10000): mmdid, n(100)

*/


cap pro drop mmdid 
program define mmdid , rclass 
	syntax , n(integer) 
	clear 
	set obs `n'
	gen id = _n
	gen treat = runiform() >= 0.5
	gen dd = cond(treat==1 , runiform(),0)
	replace dd = 1 if dd >= 0.5
	replace dd = 2 if dd > 0 & dd < 1
	
	expand 10
	sort id 
	bys id: gen t = _n 
	gen x = rnormal()
	gen ei = cond(treat==1,t - 5, -5)
	gen y = cond(ei>=0 & dd==1, 5+ t +x + rnormal() , cond(ei >=0 & dd==2 , 10 + t +x + rnormal(), t+x+rnormal())) 
	gen n = 1
	egen nj = sum(n), by(ei)
	egen njj = sum(n), by(ei dd)


	 xtmixed y i.treat##c.x i.t || ei: || dd:,

	 predict u1 u2 , reffect
	 
	 gen u11 = u1 / (exp(_b[lns1_1_1:_cons])^2 /(exp(_b[lns1_1_1:_cons]) ^2 + exp(_b[lnsig_e:_cons]) ^2 /nj))
	 gen u12 = u2 / (exp(_b[lns2_1_1:_cons])^2 /(exp(_b[lns2_1_1:_cons]) ^2 + exp(_b[lnsig_e:_cons]) ^2 /njj))
	 gen u = u11 + u12
	
	 
	 egen tag = tag(ei dd)

	 gen pp=.

	 sum u if ei>-5 & ei <0 & dd==1 & tag==1
	 replace pp = r(mean) if dd==1
	 sum u if ei>-5 & ei <0 & dd==2 & tag==1
	 replace pp = r(mean) if dd==2
	 
	 gen ate = u - pp
		

	 egen att = mean(ate)  , by(ei) 
	 egen tag1 = tag(ei)
	 mkmat att if tag1==1 & att !=.
	 
	 forv i = 5/7 {
	 	scalar att`i' = att[`i',1]
		ret scalar att`i' = att`i'
	 }
	 
	gen ttt = cond(treat==1,5,0)
 
	
	csdid y , i(id) time(t) gvar(ttt)
	estat event
	
	ret scalar cs5 = _b[g5:t_4_5]
	ret scalar cs6 = _b[g5:t_4_6]
	ret scalar cs7 = _b[g5:t_4_7]
	end
	
