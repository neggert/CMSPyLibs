from ROOT import TLorentzVector, TVector2, TVector3
from math import *
from copy import copy

class VarCalculator:
	def __init__(self):
		self.jet1=[]
		self.jet2=[]
		self.lep1=[]
		self.lep2=[]
		self.met=[]
		self.particlesSet=False
		self.mc=0.0
		self.comE=7000.0
		self.totalE=0.0

	def reset(self):
		self.jet1=[]
		self.jet2=[]
		self.lep1=[]
		self.lep2=[]
		self.met=[]
		self.particlesSet=False
		self.mc=0.0
		self.comE=7000.0
		self.totalE=0.0

	def setParticles(self, jet1in, jet2in, lep1in, lep2in, metin):
		self.jet1=TLorentzVector(jet1in.px(),jet1in.py(),jet1in.pz(),jet1in.energy())
		self.jet2=TLorentzVector(jet2in.px(),jet2in.py(),jet2in.pz(),jet2in.energy())
		self.lep1=TLorentzVector(lep1in.px(),lep1in.py(),lep1in.pz(),lep1in.energy())
		self.lep2=TLorentzVector(lep2in.px(),lep2in.py(),lep2in.pz(),lep2in.energy())
		self.met=TLorentzVector(metin.px(),metin.py(),metin.pz(),metin.energy())

		self.totalE=metin.sumEt()+metin.energy()

		self.particlesSet=True

	def setP4s(self, jet1in, jet2in, lep1in, lep2in, metin):
		self.jet1=jet1in
		self.jet2=jet2in
		self.lep1=lep1in
		self.lep2=lep2in
		self.met=metin

		self.totalE=metin.Et()+metin.E()

		self.particlesSet=True

#		print "Setting particles"
#		print self.jet1.M(), self.jet1.Px(), self.jet1.Py()
#		print self.jet2.M(), self.jet2.Px(), self.jet2.Py()
#		print self.lep1.M(), self.lep1.Px(), self.lep1.Py()
#		print self.lep2.M(), self.lep2.Px(), self.lep2.Py()
#		print self.met.M(), self.met.Px(), self.met.Py()

	def setMc(self, mcin=0.0):
		self.mc=mcin

	def setCom(self, comIn=7000.0):
		self.comE=comIn

	def w1a(self):
		return self.jet1+self.lep1

	def w2a(self):
		return self.jet2+self.lep2

	def w1b(self):
		return self.jet1+self.lep2

	def w2b(self):
		return self.jet2+self.lep1

	def child221(self):
		return self.lep1+self.lep2+self.met

	def upstream210(self):
		return -self.lep1-self.lep2-self.met

	def upstream220(self):
		return -self.jet1-self.jet2-self.lep1-self.lep2-self.met

	def upstream221(self):
		return -self.jet1-self.jet2-self.lep1-self.lep2-self.met

	def calcPerp(self,p4Vis, p4Upstream):
		scale=(p4Vis.Px()*p4Upstream.Px()+p4Vis.Py()*p4Upstream.Py())/(p4Upstream.Px()*p4Upstream.Px()+p4Upstream.Py()*p4Upstream.Py())
		perpOutput=TLorentzVector()
		perpOutput.SetXYZM(p4Vis.Px()-scale*(p4Upstream.Px()),p4Vis.Py()-scale*(p4Upstream.Py()),0.0,p4Vis.M())
		return perpOutput

	def calcPar(self,p4Vis, p4Upstream):
		scale=(p4Vis.Px()*p4Upstream.Px()+p4Vis.Py()*p4Upstream.Py())/(p4Upstream.Px()*p4Upstream.Px()+p4Upstream.Py()*p4Upstream.Py())
		parOutput=TLorentzVector()
		parOutput.SetXYZM(scale*p4Upstream.Px(),scale*p4Upstream.Py(),0.0,p4Vis.M())
		return parOutput
	
	def sMin(self) :
		if self.particlesSet :
			p4vis = self.jet1+self.jet2+self.lep1+self.lep2
			return sqrt(p4vis.M2()+p4vis.Pt()**2)+sqrt(self.met.M2()+self.met.Pt()**2)
		else :
			return 0.0

	def blInvarMass(self):
		if self.particlesSet:
			return ((self.w1a().M(),self.w2a().M()),(self.w1b().M(),self.w2b().M()))
		else:
			return ((0.0,0.0),(0.0,0.0))

	def mt2_210(self):
#		print "210"
		if self.particlesSet:
			return self.calcMt2(self.lep1, self.lep2, self.met)
		else:
			return 0.0

	def mt2Perp_210(self):
#		print "210perp"
		if self.particlesSet:
			return self.calcMt2(self.calcPerp(self.lep1,self.upstream210()),self.calcPerp(self.lep2,self.upstream210()),self.calcPerp(self.met,self.upstream210()))
		else:
			return 0.0

	def mt2Par_210(self):
#		print "210par"
		if self.particlesSet:
			return self.calcMt2(self.calcPar(self.lep1,self.upstream210()),self.calcPar(self.lep2,self.upstream210()),self.calcPar(self.met,self.upstream210()))
		else:
			return 0.0

	def mt2_220(self):
		#print "220"
		#print "jet1 ",self.jet1.Px(), self.jet1.Py(), self.jet1.Pz(), self.jet1.M(),self.jet1.E()
		#print "jet2 ",self.jet2.Px(), self.jet2.Py(), self.jet2.Pz(), self.jet2.M(),self.jet2.E()
	        #print "lep1 ",self.lep1.Px(), self.lep1.Py(), self.lep1.Pz(), self.lep1.M(),self.lep1.E()
		#print "lep2 ",self.lep2.Px(), self.lep2.Py(), self.lep2.Pz(), self.lep2.M(),self.lep2.E()
		#print "met ",self.met.Px(),self.met.Py(),self.met.Pz(),self.met.M(),self.met.E()
		#print "w1a ",self.w1a().Px(),self.w1a().Py(),self.w1a().Pz(),self.w1a().M(),self.w1a().E()
		#print "w2a ",self.w2a().Px(),self.w2a().Py(),self.w2a().Pz(),self.w2a().M(),self.w2a().E()
		#print "w1b ",self.w1b().Px(),self.w1b().Py(),self.w1b().Pz(),self.w1b().M(),self.w1b().E()
		#print "w2b ",self.w2b().Px(),self.w2b().Py(),self.w2b().Pz(),self.w2b().M(),self.w2b().E()
		
		if self.particlesSet:
			return (self.calcMt2(self.w1a(), self.w2a(), self.met),self.calcMt2(self.w1b(), self.w2b(), self.met))
		else:
			return (0.0,0.0)

	def mt2Perp_220(self):
#		print "220perp"
		if self.particlesSet:
			return (self.calcMt2(self.calcPerp(self.w1a(),self.upstream220()), self.calcPerp(self.w2a(),self.upstream220()), self.calcPerp(self.met,self.upstream220())),self.calcMt2(self.calcPerp(self.w1b(),self.upstream220()), self.calcPerp(self.w2b(),self.upstream220()), self.calcPerp(self.met,self.upstream220())))
		else:
			return (0.0,0.0)

	def mt2Par_220(self):
#		print "220par"
		if self.particlesSet:
			return (self.calcMt2(self.calcPar(self.w1a(),self.upstream220()), self.calcPar(self.w2a(),self.upstream220()), self.calcPar(self.met,self.upstream220())),self.calcMt2(self.calcPar(self.w1b(),self.upstream220()), self.calcPar(self.w2b(),self.upstream220()), self.calcPar(self.met,self.upstream220())))
		else:
			return (0.0,0.0)

	def mt2_221(self):
#		print "221"
		if self.particlesSet:
			return self.calcMt2(self.jet1, self.jet2, self.child221())
		else:
			return 0.0

	def mt2Perp_221(self):
#		print "221perp"
		if self.particlesSet:
			return self.calcMt2(self.calcPerp(self.jet1,self.upstream221()),self.calcPerp(self.jet2,self.upstream221()),self.calcPerp(self.child221(),self.upstream221()))
		else:
			return 0.0

	def mt2Par_221(self):
#		print "221par"
		if self.particlesSet:
			return self.calcMt2(self.calcPar(self.jet1,self.upstream221()),self.calcPar(self.jet2,self.upstream221()),self.calcPar(self.child221(),self.upstream221()))
		else:
			return 0.0

	def mct_210(self):
		if self.particlesSet:
			return self.calcMct(self.lep1, self.lep2)
		else:
			return 0.0

	def mctPerp_210(self):
		if self.particlesSet:
			return self.calcMct(self.calcPerp(self.lep1,self.upstream210()),self.calcPerp(self.lep2,self.upstream210()))
		else:
			return 0.0

	def mctPar_210(self):
		if self.particlesSet:
			return self.calcMct(self.calcPar(self.lep1,self.upstream210()),self.calcPar(self.lep2,self.upstream210()))
		else:
			return 0.0

	def mct_220(self):
		if self.particlesSet:
			return (self.calcMct(self.w1a(), self.w2a()),self.calcMct(self.w1b(), self.w2b()))
		else:
			return (0.0,0.0)

	def mctPerp_220(self):
		if self.particlesSet:
			return (self.calcMct(self.calcPerp(self.w1a(),self.upstream220()), self.calcPerp(self.w2a(),self.upstream220())), self.calcMct(self.calcPerp(self.w1b(),self.upstream220()), self.calcPerp(self.w2b(),self.upstream220())))
		else:
			return (0.0,0.0)

	def mctPar_220(self):
		if self.particlesSet:
			return (self.calcMct(self.calcPar(self.w1a(),self.upstream220()), self.calcPar(self.w2a(),self.upstream220())), self.calcMct(self.calcPar(self.w1b(),self.upstream220()), self.calcPar(self.w2b(),self.upstream220())))
		else:
			return (0.0,0.0)

	def mct_221(self):
		if self.particlesSet:
			return self.calcMct(self.jet1, self.jet2)
		else:
			return 0.0

	def mctPerp_221(self):
		if self.particlesSet:
			return self.calcMct(self.calcPerp(self.jet1,self.upstream221()),self.calcPerp(self.jet2,self.upstream221()))
		else:
			return 0.0

	def mctPar_221(self):
		if self.particlesSet:
			return self.calcMct(self.calcPar(self.jet1,self.upstream221()),self.calcPar(self.jet2,self.upstream221()))
		else:
			return 0.0

	def mctBoostCor_210(self):
		if self.particlesSet:
			return self.calcMctBoostCor(self.lep1, self.lep2, self.upstream210())
		else:
			return (0.0,0)

	def mctBoostCor_220(self):
		if self.particlesSet:
			return (self.calcMctBoostCor(self.w1a(),self.w2a(),self.upstream220()), self.calcMctBoostCor(self.w1b(),self.w2b(),self.upstream220()) )
		else:
			return ((0.0,0),(0.0,0))

	def mctBoostCor_221(self):
		if self.particlesSet:
			return self.calcMctBoostCor(self.jet1, self.jet2, self.upstream221())
		else:
			return (0.0,0)

	def calcMt2(self, vis1in, vis2in, childin):
#		print "Input"
#		print vis1in.M(),vis1in.Px(),vis1in.Py()
#		print vis2in.M(),vis2in.Px(),vis2in.Py()
		if vis1in.M()>vis2in.M():
			vis1=copy(vis2in)
			vis2=copy(vis1in)
		else:
			vis1=copy(vis1in)
			vis2=copy(vis2in)

		child=copy(childin)

		mag=self.mt2Sqrt(vis1.Px()*vis1.Px()+vis1.Py()*vis1.Py())
		cospart=vis1.Px()/mag
		sinpart=vis1.Py()/mag

		ma,pxa,pya=vis1.M(),vis1.Px(),vis1.Py()
		mb,pxb,pyb=vis2.M(),vis2.Px(),vis2.Py()
		mc,pxc,pyc=child.M(),child.Px(),child.Py()

		vis1.SetXYZM(mag,0.0,0.0,ma)
		vis2.SetXYZM(pxb*cospart+pyb*sinpart, pyb*cospart-pxb*sinpart, 0.0, mb)
		child.SetXYZM(pxc*cospart+pyc*sinpart, pyc*cospart-pxc*sinpart, 0.0, mc)

#		print "After rotating, new copies"
#		print vis1.M(),vis1.Px(),vis1.Py()
#		print vis2.M(),vis2.Px(),vis2.Py()

#		print "Supposed set mass of copies"
#		print ma,mb

#		print "New particles"
#		print vis1.M(),vis1.Px(),vis1.Py()
#		print vis2.M(),vis2.Px(),vis2.Py()
#		print self.mc,child.Px(),child.Py()

		outputmt2=0.0
		solved=False

		vis1M2=vis1.M()*vis1.M()
		vis2M2=vis2.M()*vis2.M()
		mc2=self.mc*self.mc
		vis1Px2=vis1.Px()*vis1.Px()
		vis2Px2=vis2.Px()*vis2.Px()
		childPx2=child.Px()*child.Px()
		vis2Py2=vis2.Py()*vis2.Py()
		childPy2=child.Py()*child.Py()
		vis1Pt2=vis1Px2
		vis2Pt2=vis2Px2+vis2Py2
		childPt2=childPx2+childPy2
		vis1Et2=vis1M2+vis1Pt2
		vis2Et2=vis2M2+vis2Pt2
		childEt2=mc2+childPt2
		vis1Et=self.mt2Sqrt(vis1Et2)
		vis2Et=self.mt2Sqrt(vis2Et2)

		Mtmin=0.0
		Mtmax=0.0

		if (not(vis1.M()<=0.0 or vis2.M()<=0.0)):
#			print "This should not give division errors"
#			print vis1.M(),vis2.M()
			xlmin=vis1.Px()*self.mc/vis1.M()
			xrmin=vis2.Px()*self.mc/vis2.M()
			yrmin=vis2.Py()*self.mc/vis2.M()

			altxlmin=child.Px()-xlmin
			altxrmin=child.Px()-xrmin
			altyrmin=child.Py()-yrmin

			Mtlmin=vis1.M()+self.mc
			Mtrmin=vis2.M()+self.mc

			Mtratlmin=self.mt2Sqrt(vis2M2+mc2+2.0*(vis2Et*self.mt2Sqrt(mc2+altxlmin*altxlmin+childPy2)-vis2.Px()*altxlmin-vis2.Py()*child.Py()))
			Mtlatrmin=self.mt2Sqrt(vis1M2+mc2+2.0*(vis1Et*self.mt2Sqrt(mc2+altxrmin*altxrmin+altyrmin*altyrmin)-vis1.Px()*altxrmin))

			if (Mtlmin>=Mtratlmin):
				solved=True
				outputmt2=Mtlmin
#				print "New set at point 1"
			elif (Mtrmin>=Mtlatrmin):
				solved=True
				outputmt2=Mtrmin
#				print "New set at point 2"
			else:
				if (Mtlmin>Mtrmin):
					Mtmin=Mtlmin
				else:
					Mtmin=Mtrmin
				
				if (Mtlatrmin<Mtratlmin):
					Mtmax=Mtlatrmin
				else:
					Mtmax=Mtratlmin

			backupmid=self.mt2Sqrt(mc2+(child.Px()*child.Px()+child.Py()*child.Py())*0.25)
			backup1=self.mt2Sqrt(vis1M2+mc2+2.0*(vis1Et*backupmid-0.5*vis1.Px()*child.Px()))
			backup2=self.mt2Sqrt(vis2M2+mc2+2.0*(vis2Et*backupmid-0.5*(vis2.Px()*child.Px()+vis2.Py()*child.Py())))
			if (backup1>backup2):
				if (backup1<Mtmax):
					Mtmax=backup1
			else:
				if (backup2<Mtmax):
					Mtmax=backup2

		elif (not(vis2.M()<=0.0)):
			xrmin=vis2.Px()*self.mc/vis2.M()
			yrmin=vis2.Py()*self.mc/vis2.M()

			altxrmin=child.Px()-xrmin
			altyrmin=child.Py()-yrmin

			Mtlmin=self.mc
			Mtrmin=vis2.M()+self.mc

			Mtlatrmin=self.mt2Sqrt(vis1M2+mc2+2.0*(vis1Et*self.mt2Sqrt(mc2+altxrmin*altxrmin+altyrmin*altyrmin)-vis1.Px()*altxrmin))

			if (Mtrmin>=Mtlatrmin):
				solved=True
				outputmt2=Mtrmin
#				print "New set at point 3"
			else:
				if (Mtlmin>Mtrmin):
					Mtmin=Mtlmin
				else:
					Mtmin=Mtrmin

				Mtmax=Mtlatrmin

			backupmid=self.mt2Sqrt(mc2+(child.Px()*child.Px()+child.Py()*child.Py())*0.25)
			backup1=self.mt2Sqrt(vis1M2+mc2+2.0*(vis1Et*backupmid-0.5*vis1.Px()*child.Px()))
			backup2=self.mt2Sqrt(vis2M2+mc2+2.0*(vis2Et*backupmid-0.5*(vis2.Px()*child.Px()+vis2.Py()*child.Py())))
			if (backup1>backup2):
				if (backup1<Mtmax):
					Mtmax=backup1
			else:
				if (backup2<Mtmax):
					Mtmax=backup2

		else:
			Mtmin=self.mc
			trialmid=self.mt2Sqrt(mc2+0.25*childPt2)
			trial1=self.mt2Sqrt(mc2+2.0*(vis1Et*trialmid-vis1.Px()*child.Px()*0.5))
			trial2=self.mt2Sqrt(mc2+2.0*(vis2Et*trialmid-0.5*(vis2.Px()*child.Px()+vis2.Py()*child.Py())))
			if (trial1>trial2):
				Mtmax=trial1
			else:
				Mtmax=trial2

#		print "New bounds",Mtmin,Mtmax

		if (not solved):
			solved=True
			outputmt2=Mtmin+(Mtmax-Mtmin)*0.5

			C1=1.0/vis1Et2

			A2=vis2M2+vis2Py2
			B2=-2.0*vis2.Px()*vis2.Py()
			C2=1.0/(vis2M2+vis2Px2)

			preF1=vis1Et2*mc2

			preD2=-2.0*child.Px()*A2-B2*child.Py()
			preE2=-2.0*child.Py()/C2-B2*child.Px()
			preF2=vis2Et2*childEt2-childPx2*vis2Px2-childPy2*vis2Py2+B2*child.Px()*child.Py()

			G=B2*0.5*C2
			J1=-vis1M2*C1
			J2=(B2*B2*0.25*C2-A2)*C2

			alpha=G*G-J1-J2
			p0_4=alpha*alpha-4.0*J1*J2

			while (outputmt2>Mtmin and outputmt2<Mtmax):
				q1=mc2+vis1M2-outputmt2*outputmt2
				D1=q1*vis1.Px()
				F1=preF1-q1*q1*0.25

				q2=outputmt2*outputmt2-mc2-vis2M2
				D2=preD2+q2*vis2.Px()
				E2=preE2+q2*vis2.Py()
				F2=preF2-q2*(q2*0.25+vis2.Px()*child.Px()+vis2.Py()*child.Py())

				H=E2*0.5*C2

				K1=-D1*C1
				L1=-F1*C1

				K2=(B2*E2*0.5*C2-D2)*C2
				L2=(E2*E2*0.25*C2-F2)*C2

				beta=2.0*G*H-K1-K2
				gamma=H*H-L1-L2

				p0_4nonzero=p0_4
				if fabs(p0_4)<1e-9:
					p0_4nonzero=1e-9

				p0_3=(2.0*alpha*beta-4.0*(J1*K2+J2*K1))/p0_4nonzero
				p0_2=(2.0*alpha*gamma+beta*beta-4.0*(J1*L2+J2*L1+K1*K2))/p0_4nonzero
				p0_1=(2.0*beta*gamma-4.0*(K1*L2+K2*L1))/p0_4nonzero
				p0_0=(gamma*gamma-4.0*L1*L2)/p0_4nonzero

				p2_2=0.1875*p0_3*p0_3-p0_2*0.5
				p2_1=p0_3*p0_2*0.125-0.75*p0_1
				p2_0=p0_3*p0_1*0.0625-p0_0

				p2_2nonzero=p2_2
				if fabs(p2_2)<1e-9:
					p2_2nonzero=1e-9

				p3_1=(4.0*p2_0+3.0*p0_3*p2_1)/p2_2nonzero-4.0*p2_1*p2_1/(p2_2nonzero*p2_2nonzero)-2.0*p0_2
				p3_0=3.0*p0_3*p2_0/p2_2nonzero-4.0*p2_1*p2_0/(p2_2nonzero*p2_2nonzero)-p0_1

				p4_0=0
				if fabs(p3_1)<1e-9:
					if p2_2>0:
						p4_0=-1
					elif p2_2==0:
						p4_0=0
					else:
						p4_0=1
				else:
					p4_0=p2_1*p3_0/p3_1-p2_2*p3_0*p3_0/(p3_1*p3_1)-p2_0

				negroots=1
				posroots=0

				if ((p0_4<0.0 and p2_2<0.0) or (p0_4>0.0 and p2_2>0.0)):
					negroots+=1
				elif((p0_4<0.0 and p2_2>0.0) or (p0_4>0.0 and p2_2<0.0)):
					posroots+=1

				if ((p2_2<0.0 and p3_1<0.0) or (p2_2>0.0 and p3_1>0.0)):
					negroots+=1
				elif((p2_2<0.0 and p3_1>0.0) or (p2_2>0.0 and p3_1<0.0)):
					posroots+=1

				if ((p3_1<0.0 and p4_0<0.0) or (p3_1>0.0 and p4_0>0.0)):
					negroots+=1
				elif((p3_1<0.0 and p4_0>0.0) or (p3_1>0.0 and p4_0<0.0)):
					posroots+=1

				if (posroots==negroots):
					Mtmin=outputmt2
				else:
					Mtmax=outputmt2

				outputmt2=Mtmin+(Mtmax-Mtmin)*0.5
#		print "New set at point 4"
		return outputmt2

	def calcMct(self, p4Vis1, p4Vis2):
		Et1 = self.mt2Sqrt(p4Vis1.Perp2()+p4Vis1.M2())
		Et2 = self.mt2Sqrt(p4Vis2.Perp2()+p4Vis2.M2())

		try :
			mct = self.mt2Sqrt(p4Vis1.M2()+p4Vis2.M2()+2.0*(Et1*Et2+p4Vis1.Px()*p4Vis2.Px()+p4Vis1.Py()*p4Vis2.Py()))
		except ValueError :
			return 0.
		if isnan(mct) :
			return 0.
		return mct

	def calcMctBoostCor(self, p4Vis1, p4Vis2, p4Upstream):    # for totalE, suggest met.sumet()+met.energy()
		if (self.totalE > self.comE) :
			print "totalE is greater than comEnergy"

		# rotate such that upstream momentum is in x -direction
		p4Vis1Rot = copy(p4Vis1)
		p4Vis2Rot = copy(p4Vis2)
		p4UpstreamRot = copy(p4Upstream)

		phi = atan2(p4Upstream.Py(),p4Upstream.Px())
		p4UpstreamRot.RotateZ(-phi)
		p4Vis1Rot.RotateZ(-phi)
		p4Vis2Rot.RotateZ(-phi)

		try :
			Ax = p4Vis1Rot.Px()*self.mt2Sqrt(p4Vis2Rot.Py()**2+p4Vis2Rot.M2())+p4Vis2Rot.Px()*self.mt2Sqrt(p4Vis1Rot.Py()**2+p4Vis1Rot.M2())
		except ValueError :
			return (0., -1)

		# boost the visible particles with the smallest possible boost
		p4Vis1Lo = copy(p4Vis1Rot)
		p4Vis1Lo.Boost(p4UpstreamRot.Px()/self.comE,0.,0.)
		p4Vis2Lo = copy(p4Vis2Rot)
		p4Vis2Lo.Boost(p4UpstreamRot.Px()/self.comE,0.,0.)

		if (Ax > 0.0 ) :
			return (self.calcMct(p4Vis1Lo, p4Vis2Lo), 0)

		try :
			AxLo = p4Vis1Lo.Px()*self.mt2Sqrt(p4Vis2Lo.Py()**2+p4Vis2Lo.M2())+p4Vis2Lo.Px()*self.mt2Sqrt(p4Vis1Lo.Py()**2+p4Vis1Lo.M2())
		except ValueError :
			return (0., -1.)

		if (AxLo > 0) :
			return (self.calcMct(p4Vis1Lo, p4Vis2Lo), 1)


		# boost the particles with the largest possible boost
		if self.totalE == 0. :
			return (0., -1)
		p4Vis1Hi = copy(p4Vis1Rot)
		p4Vis1Hi.Boost(p4UpstreamRot.Px()/self.totalE,0.,0.)
		p4Vis2Hi = copy(p4Vis2Rot)
		p4Vis2Hi.Boost(p4UpstreamRot.Px()/self.totalE,0.,0.)

		try : 
			AxHi = p4Vis1Hi.Px()*self.mt2Sqrt(p4Vis2Hi.Py()**2+p4Vis2Hi.M2())+p4Vis2Hi.Px()*self.mt2Sqrt(p4Vis1Hi.Py()**2+p4Vis1Hi.M2())
		except ValueError :
			return (0., -1)

		if (AxHi < 0.0 ) :
			return (self.calcMct(p4Vis1Hi, p4Vis2Hi), 2)
		else :
			return (self.calcMct(self.calcPerp(p4Vis1,p4Upstream), self.calcPerp(p4Vis2,p4Upstream)), 3)


	def mt2Sqrt(self,x):
		if isnan(x):
			return 0.0
		elif (x<=0.0):
			return 0.0
		elif (isinf(x)):
			return 1e99999999999999999999999999999999
		else:
			prev2root=-1.0
			prevroot=-1.0
			root=1.0
			while((root!=prevroot) and (root!=prev2root)):
				prev2root=prevroot
				prevroot=root
				root=(root+x/root)*0.5

			return root
