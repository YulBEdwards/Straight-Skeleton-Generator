import matplotlib.pyplot as plt
import math
import copy

# Line segment class
class Segment:
    def __init__(self, sid, tid, nextid, previd, origid, p1, p2, a, b, c):
        self.sid    = sid  # Unique identifier
        self.tid    = tid
        self.nextid = nextid
        self.previd = previd
        self.origid = origid
        self.p1     = p1   # (x1, y1)
        self.p2     = p2   # (x2, y2)
        self.a      = a    # line is a*x+b*y+c
        self.b      = b
        self.c      = c

class Line:
    def __init__(self, a, b, c, p1, p2):
        self.p1  = p1
        self.p2  = p2
        self.a   = a    # line is a*x+b*y+c
        self.b   = b
        self.c   = c

class Ray:
    def __init__(self, a, b, c, p):
        self.p   = p    # p is the origin of the ray
        self.a   = a    # line is a*x+b*y+c
        self.b   = b
        self.c   = c
    
class Track:
    def __init__(self, tid, sid, leaf, parting, t, speed, \
                 bisect, target, clipid, caped, niid):
        self.tid      = tid     # Track identifier
        self.sid      = sid     # Segment identifier
        self.leaf     = leaf    # Side of track for continuation
        self.parting  = parting
        self.t        = t       # Time at ending intersection
        self.speed    = speed
        self.bisect   = bisect  # Bisector definition for track (p1/p2,a,b,c)
        self.target   = target  # Side id of target        
        self.clipid   = clipid  # Track that clips this track
        self.caped    = caped   # caped if intersected or through
        self.niid     = niid    # track id of next island track (or None)
        
# Track intersect object (track intersection records)
class TI:
    def __init__(self, idsml, idlrg, tsml, tlrg, xi, yi):
        self.idsml = idsml
        self.idlrg = idlrg
        self.tsml  = tsml
        self.tlrg  = tlrg
        self.xi    = xi
        self.yi    = yi

class Region:       
    def __init__(self, rid, startid, endid, starttid, endtid):
        self.rid       = rid
        self.startid   = startid  # segment indices
        self.endid     = endid
        self.starttid  = starttid # track indices       
        self.endtid    = endtid
        
class Skeleton:
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
        
# Skeleton element (Used in skelton for 3 or more sides)
class SE:
    def __init__(self, pedge, nedge):
        self.pedge = pedge
        self.nedge = nedge
#___________________________________________________________
        
def generateSegments():
    global title
        # Berg
    title="Berg"
    print("Berg")
    s = [[858, -341], [650, -427],[612, -104],[490, -177],
         [368, -98],  [203, -204],[309, -334],[272, -483],
         [152, -384], [108, -663],[245, -841],[389, -752],
         [614, -1011],[519, -599],[738, -667],[858, -341]]

#___________________________________________________________
    
    Segments = []         
    numSegments = len(s)-1
    clockwise = 0
    for i in range(numSegments+2):
        if i>numSegments-1:
            ix = i-numSegments
        else:
            ix=i
        p1 = s[ix]
        p2 = s[ix+1]
        if i<numSegments:
            cross = p1[0]*p2[1]-p1[1]*p2[0]
            clockwise = clockwise+cross
        a  = p2[1]-p1[1]
        b  = p1[0]-p2[0]
        xm = 0.5*(p1[0]+p2[0])
        ym = 0.5*(p1[1]+p2[1])
        c  = -a*xm-b*ym
        L  = math.sqrt(a**2+b**2)
        if L>0:
            a = a/L; b = b/L; c = c/L
        n = (i+1)%numSegments
        p = (i-1)%numSegments
        Segments.append(Segment(sid=i, tid=None, nextid=n, previd=p,
                                origid=i, p1=p1, p2=p2, a=a, b=b, c=c))
        
    if clockwise<0:
        print("closes clockwise")
        Print("Convention requires counter clockwise closure")
    else:
        print("closes counter clockwise")
        
        
    return Segments, numSegments

def generateTracks(Segments):
    corners = []
    tracks  = []
    num_Corners = numSegments   
    clockwise = 0
    for i in range(num_Corners):
        s1 = Segments[i]
        s2 = Segments[i+1]
        turning = s1.a*s2.b-s1.b*s2.a
        reflex = turning<0
        if reflex:
                # determine bisector
            a = s2.a-s1.a; b = s2.b-s1.b; c = s2.c-s1.c
            L  = math.sqrt(a**2+b**2)
            if L>0:
                a = a/L; b = b/L; c = c/L
            bisect = Segment(sid=i, tid=None, nextid=None, previd=None,
                             origid=None, p1=s1.p2, p2=[0,0], a=a, b=b, c=c)   
            #bisect = Segment(p1=s1.p2, p2=[0,0], a=a, b=b, c=c)
            tid = len(tracks)
            Segments[i].tid = len(tracks)      
            cos_turning = s1.a*s2.a+s1.b*s2.b
            speed = math.sin(0.5*(math.acos(cos_turning)+math.pi))            
            tracks.append(Track(tid = tid, sid=i, leaf = None, t = 0, \
                                parting=True, speed=speed, \
                                bisect=bisect, target=0, clipid=None, \
                                caped=False, niid=None))
            
    return corners, tracks, Segments

# find target of a ray starting at vertex corner
def intersectionSegmentTrack(segment,track):
    if track.sid==numSegments-1:
        ix = 0
    else:
        ix = track.sid+1
    if segment.sid == track.sid or segment.sid == ix:        
        return False, None, None, None
    x1 = segment.p1[0]; y1 = segment.p1[1]
    x2 = segment.p2[0]; y2 = segment.p2[1]
    x0 = track.bisect.p1[0]; y0 = track.bisect.p1[1]
    dx = track.bisect.b; dy = -track.bisect.a
    delta = dx*(y2-y1)-dy*(x2-x1)
    t = ((x1-x0)*(y2-y1)-(y1-y0)*(x2-x1))/delta
    u = -(dx*(y1-y0)-dy*(x1-x0))/delta
    intersect  = u>=0 and u<=1 and t>=0
    xi = x0+t*dx; yi = y0+t*dy
    target = segment.sid    
    # get the closest segment to the base point
    if intersect and track.bisect.p2 != [0,0]:
        d1 = (track.bisect.p1[0]-track.bisect.p2[0])**2+ \
             (track.bisect.p1[1]-track.bisect.p2[1])**2
        d2 = (track.bisect.p1[0]-xi)**2+ \
             (track.bisect.p1[1]-yi)**2
        if d1<d2:         
            return False, None, None, None
    return intersect,xi,yi,target

def intersectionTrackTrack(track1,track2):
    x1 = track1.bisect.p1[0]; y1 = track1.bisect.p1[1]
    x2 = track1.bisect.p2[0]; y2 = track1.bisect.p2[1]
    x3 = track2.bisect.p1[0]; y3 = track2.bisect.p1[1]
    x4 = track2.bisect.p2[0]; y4 = track2.bisect.p2[1]
    t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/  \
        ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
    u = -((x1-x2)*(y1-y3)-(y1-y2)*(x1-x3))/  \
        ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
    t = round(t,4); u = round(u,4)
    intersect  = u<=1 and u>=0 and t<=1 and t>=0
    xi = x1+t*(x2-x1); yi = y1+t*(y2-y1)
    return intersect,xi,yi

def distanceToEdge(p,segment):
    return p[0]*segment.a+p[1]*segment.b+segment.c

    # Generate bisector of two segments
def generateBisector(seg1,seg2):
    a = seg2.a-seg1.a; b = seg2.b-seg1.b; c = seg2.c-seg1.c
    L  = math.sqrt(a**2+b**2)
    if L>0:
        a = a/L; b = b/L; c = c/L
            # for consecutive lines, p is the bisector origin
    return Line(a=a, b=b, c=c, p1=None, p2 = None)

    # Generate bisector of two segments
def generateRay(seg1,seg2):
    a = seg2.a-seg1.a; b = seg2.b-seg1.b; c = seg2.c-seg1.c
    L  = math.sqrt(a**2+b**2)
    if L>0:
        a = a/L; b = b/L; c = c/L
            # for consecutive lines, p is the bisector origin
    return Ray(a=a, b=b, c=c, p=seg1.p2)

    # Intersection of two tracks
def intersectionTrackBisector(track,bisect):
    x1 = track.bisect.p1[0]; y1 = track.bisect.p1[1]
    x2 = track.bisect.p2[0]; y2 = track.bisect.p2[1]
    x0 = bisect.p[0]; y0 = bisect.p[1]
    dx = bisect.b; dy = -bisect.a
    delta = dx*(y2-y1)-dy*(x2-x1)
    t = ((x1-x0)*(y2-y1)-(y1-y0)*(x2-x1))/delta
    u = -(dx*(y1-y0)-dy*(x1-x0))/delta
    t = round(t,4); u = round(u,4)
    intersect  = u>=0 and u<=1 and t<=0
    xi = x1+u*(x2-x1); yi = y1+u*(y2-y1)
    return intersect,xi,yi

    # Intersection of bisectors of two pairs of segments
def intersectionBisectBisect(seg1A,seg1B,seg2A,seg2B):
    LD1 = generateRay(Segments[seg1A],Segments[seg1B])
    LD2 = generateRay(Segments[seg2A],Segments[seg2B])
    delta = LD1.a*LD2.b-LD2.a*LD1.b
    x = (LD1.b*LD2.c-LD2.b*LD1.c)/delta
    y = (LD2.a*LD1.c-LD1.a*LD2.c)/delta
    return x,y

    # Find an island among successors of clipped track
def findIsland(region,pass2):
    startid = region.endtid # startid and endtid are track ids
    curid = startid
    if Tracks[curid].clipid==None:
        return False
    for i in range(len(Tracks)):
        if pass2:
            Tracks[curid].niid = Tracks[curid].clipid
        nextid=Tracks[curid].clipid    
        if startid==nextid:
            return True
        curid = Tracks[curid].clipid
        if curid==None:
            return False
    return False

def findTrueTarget():
    for track in Tracks:        
        if track.clipid == None: # If through track (clipid=None), cap it
            x,y = intersectionBisectBisect(track.sid,(track.sid+1)&numSegments, \
                                           track.sid,track.target)
            track.bisect.p2 = [x,y]            
            track.caped = True
            continue
        
        if not track.caped:
            clippingTrack = Tracks[track.clipid]
            t1= clippingTrack.sid   # the next 4 are the sides 
            t2 = (t1+1)%numSegments
            s1 = track.sid
            s2 = (s1+1)%numSegments
            d = distanceToEdge(track.bisect.p1,clippingTrack.bisect)
            niid = Tracks[clippingTrack.tid].niid
            if track.niid!=None:    # If not None, part of an island
                 # Examine all for & aft segments at track bases
                niid = Tracks[niid].tid
                u1 = Tracks[niid].sid
                u2 = (u1+1)%numSegments
                if track.tid==5:
                    niid = 1
                if track.tid==1:
                    niid = 4
                if track.tid==4:
                    niid = 1
                u1 = Tracks[niid].sid
                u2 = (u1+1)%numSegments
                x1,y1 = intersectionBisectBisect(s1,s2,s1,u1)          
                x2,y2 = intersectionBisectBisect(s1,s2,s1,u2)
                ax.plot([x1],[y1], 'o', color ="red", ms=3)
                ax.plot([x2],[y2], 'o', color ="black", ms=3)
            if d>=0:
                x,y = intersectionBisectBisect(s1,s2,s1,t1)
                track.target = clippingTrack.sid
            else:
                x,y = intersectionBisectBisect(s1,s2,s1,t2)
                track.target = (clippingTrack.sid+1)%numSegments
            track.bisect.p2 = [x,y]
            track.caped = True
            track.t = 0
    return

def generateRegions():
    Regions = []
    for track in Tracks:
        if track==Tracks[-1]:
            nextTrack = Tracks[0]
        else:
            nextTrack = Tracks[Tracks.index(track)+1]
        nextside = track.sid+1
        if nextside == numSegments:
            nextside = 0
        Regions.append(Region(rid=track.tid,  \
            startid=track.sid+1, endid=nextTrack.sid, \
            starttid=track.tid, endtid=nextTrack.tid))
    return Regions

 # find the starting side for a polygon with 3 or more sides
     # if istart<iend or istart>iend
     #     the region is 3 or more sides between tracks
     # If istart=iend the region is a convex polygon
def findFirst(region):
    istart = region.startid
    iend   = region.endid
    if iend<istart:
        iend = iend+numSegments    
    if iend-istart==2:
        i0 = Segments[istart].nextid
        return i0
    d0=0
    side = istart
    while True:
        mside=side%numSegments
        pside = Segments[side].previd
        #if mside==region.startid:
        #    pside=region.endid
        nside = Segments[side].nextid
        #if mside==region.endid:
        #    nside=region.startid
        xi,yi = intersectionBisectBisect(pside,mside,mside,nside)
        d = abs(distanceToEdge([xi,yi],Segments[mside]))
        if d<d0 or d0==0:
            d0=d; i0=mside; xi0=xi; yi0=yi
        side = Segments[side].nextid
        if (side>iend and istart!=iend) or side==istart:
            break
    return i0

def generateSkeletonSegment(p1,p2):
        x1 = p1[0]; y1 = p1[1]
        x2 = p2[0]; y2 = p2[1]
        Skeletons.append(Skeleton(p1=[x1,y1],p2=[x2,y2]))
        ax.plot([x1,x2],[y1,y2],
            color = "green", linewidth=0.75 )
        ax.plot([x2],[y2], 'o', color ="green", ms=1.75)

# Spine march is for three adjacent convex sides
def spineMarch(n,region):  # n is the segment of the base of first
        # Order of the edges in seA and seB:
        #  pedge seA nedge ... pedge seB nedge
    x,y=intersectionBisectBisect((n-1)%numSegments,n,n,(n+1)%numSegments)
    generateSkeletonSegment(Segments[n].p1,[x,y])
    generateSkeletonSegment(Segments[n].p2,[x,y])
    
    seA=SE(pedge=(n-2)%numSegments,nedge=(n-1)%numSegments)
    seB=SE(pedge=(n+1)%numSegments,nedge=(n+2)%numSegments)
    trail=[seA.nedge,seB.pedge]
    while True:
        xia,yia = intersectionBisectBisect(seA.pedge, seA.nedge,
                                           trail[0],trail[1])
        xib,yib = intersectionBisectBisect(seB.pedge, seB.nedge,
                                           trail[0],trail[1])      
        da = (xia-x)**2+(yia-y)**2
        db = (xib-x)**2+(yib-y)**2        
        if round(da,6)<=round(db,6):
            # Doing intersection wih A side
            Aside = True
            generateSkeletonSegment([x,y],[xia,yia])
            generateSkeletonSegment([xia,yia],Segments[trail[0]].p1)
            
            x=xia; y=yia
                # Update the edges
            seA.nedge=seA.pedge
            seA.pedge=Segments[seA.pedge].previd
                
        else:
            Aside = False
            generateSkeletonSegment([x,y],[xib,yib])
            generateSkeletonSegment([xib,yib],Segments[trail[1]].p2)
                # Set up the back coord for next step
            x=xib
            y=yib
            
                     
        c1 = Segments[seA.pedge].tid!=None
        c2 = Segments[seB.pedge].tid!=None
        
        if c1 and Tracks[Segments[seA.pedge].tid].target==seB.nedge:
            track = Tracks[Segments[seA.pedge].tid]
            c2 = False
            break
        if c2 and Tracks[Segments[seB.pedge].tid].target==seA.nedge:
            track = Tracks[Segments[seB.pedge].tid]
            c1 = False
            break
        
        if not Aside:
            # Update BSide data
            seB.pedge=seB.nedge
            if seB.nedge==region.endid:
                seB.nedge=region.startid
            else:
                seB.nedge=Segments[seB.nedge].nextid
            trail[1]=seB.pedge
           
    # Fix the links
    if c2:  # Reflex corner beyond target
        Ends.append([[seA.nedge,seB.nedge],Segments[seB.pedge].tid,[x,y]])
        m = Segments[seA.nedge].nextid
        Segments[m].previd = seB.pedge
        Segments[seB.pedge].nextid = m
        Segments[seA.nedge].nextid = seB.nedge
        Segments[seB.nedge].previd = seA.nedge
        startSeg = seB.nedge
    if c1:  # Reflex corner before target

        Segments[seA.pedge].nextid = seB.nedge
        Segments[seB.nedge].previd = seA.pedge
        Segments[seA.nedge].previd = seB.pedge
        Segments[seB.pedge].nextid = seA.nedge
        
        xia= Tracks[region.starttid].bisect.p2[0]
        yia= Tracks[region.starttid].bisect.p2[1]
        Ends.append([[seA.pedge,seB.nedge],Segments[seA.pedge].tid,[xia,yia]])
        generateSkeletonSegment([x,y],[xia,yia])
        
        x= Tracks[region.starttid].bisect.p1[0]
        y= Tracks[region.starttid].bisect.p1[1]        
        generateSkeletonSegment([x,y],[xia,yia])
        startSeg = seA.pedge
        
    return startSeg

def spikeRemoval(n): # n is the beginning side index of the spike
    v1 = Tracks[Segments[(n-1)%numSegments].tid]
    v2 = Tracks[Segments[(n+1)%numSegments].tid]        
    bs = generateRay(Segments[n],Segments[(n+1)%numSegments])
    cross1,xi1,yi1 = intersectionTrackBisector(v1,bs)
    cross2,xi2,yi2 = intersectionTrackBisector(v2,bs)
    if cross1 and not cross2:
        # Intersects the lower numbered side
        p=[xi1,yi1]
        m = Segments[n].previd        
        Ends.append([[Segments[n].previd,Segments[n].nextid],Segments[m].tid,p])
        generateSkeletonSegment(v1.bisect.p1,p)        
        Skeletons.append(Skeleton(p2=p,p1=bs.p))
        v1.bisect.p2 = p

        Segments[m].nextid = Segments[n].nextid
        Segments[Segments[n].nextid].previd = m
        Segments[n].nextid=n
        Segments[n].previd=n
        v1.caped = True
        startSeg = Segments[m].nextid
        
    if not cross1 and cross2:
        # Intersect the higher numbered side      
        p=[xi2,yi2]
        m = Segments[n].nextid        
        Ends.append([[n,Segments[Segments[n].nextid].nextid],Segments[m].tid,p])        
        generateSkeletonSegment(v2.bisect.p1,p)                
        Skeletons.append(Skeleton(p2=p,p1=bs.p))
        v2.bisect.p2 = p
        Segments[Segments[m].nextid].previd = n
        Segments[n].nextid = Segments[m].nextid
        Segments[m].nextid = m
        Segments[m].previd = m
        v2.caped = True
        startSeg = n
    
    # if both edges are associated with tracks - quit   
    if cross1 or cross2:
        ax.plot([Skeletons[-1].p1[0], Skeletons[-1].p2[0]], \
            [Skeletons[-1].p1[1], Skeletons[-1].p2[1]],  \
            color = "blue", linewidth=0.75 )
        ax.plot([p[0]],[p[1]], 'o', color ="red", ms=1.75)

    if not cross1 and not cross2:
        print("Neither side intersection for spike", n)
    
    return startSeg

def findEnd(trialEnd):
    for end in Ends:
        if trialEnd[0]==end[0][0] and trialEnd[1]==end[0][1]:
            return end[2]
    return None

def checkParting(trialPair):
    for part in Parts:
        if trialPair[0]==part[0] and trialPair[1]==part[1]:
            return part[2]
    return None

def processConvexAreas():
    s = startSeg
    # This now gets only one reflex angle
    for i in range(10):
        s1 = Segments[s]
        s2 = Segments[Segments[s].nextid]
        turning = s1.a*s2.b-s1.b*s2.a
        reflex = turning<0
            # Relink at reflex corner
        if reflex:
            nextSeg = Segments[s].nextid
            trk = Tracks[Segments[s].tid]
            last = Tracks[-1]
            copyBisect = generateBisector(Segments[s],Segments[nextSeg])
            copyBisect.p1 = trk.bisect.p1
            copyBisect.p2 = trk.bisect.p2
                # note that we retain the original track and side ids
            Tracks.append(Track(tid=trk.tid,sid=trk.sid, \
                            leaf=trk.leaf, parting=trk.parting, t=trk.t, \
                            speed=trk.speed, bisect = copyBisect, \
                            target=trk.target, clipid=trk.clipid,\
                            caped=trk.caped, niid=trk.niid))
                # Make a new target for the new track
            Tseg = Segments[trk.target] # Target segment of trk

            Segments.append(Segment(len(Segments),Tseg.tid,Tseg.nextid,Tseg.previd, \
                                    Tseg.sid,Tseg.p1,Tseg.p2,Tseg.a,Tseg.b,Tseg.c))
            trk2  = Tracks[-1]
            Tseg2 = Segments[-1]
            # Fix the ids and the links
                # fix upper region
            Tseg.nextid = Segments[trk.sid].nextid
            Segments[Segments[trk.sid].nextid].previd = Tseg.sid
                # fix the lower region
            Tseg2.previd = trk.sid
            Segments[trk.sid].nextid = Tseg2.sid
            
            trk2.tid = len(Tracks)
            trk2.target = len(Segments)-1
            trk.tid = None  # The upper region division starts at target
            r1Start = Segments.index(Tseg)
            r2Start = Segments.index(Tseg2)
            Parts.append([trk.target,Segments[trk.target].nextid,trk.bisect.p2]) #Upper
            Parts.append([trk.sid,Segments[trk2.target].origid,trk.bisect.p2])  #Lower
                # trk/seg are the old copies - trk2 and seg2 are the new copies
        s = Segments[s].nextid
        if s == startSeg:
                break

    for j in range(2):
        if j==0:
            start=r1Start
        else:
            start=r2Start
        T = Region(rid=0,startid=start,endid = start,
            starttid=5,endtid=5)
        i0 = findFirst(T)
        pSeg = Segments[i0].previd; nSeg = Segments[i0].nextid
        pair1 = [pSeg,i0]; pair2 = [i0,Segments[Segments[i0].nextid].origid]
        x,y = intersectionBisectBisect(pSeg,i0,i0,nSeg)  
        p1 = findEnd(pair1)
        generateSkeletonSegment(p1,[x,y])
        p2 = findEnd(pair2)
        if p2 == None:
            ray = generateRay(Segments[pair2[0]],Segments[pair2[1]])
            p2 = ray.p
            p2 = checkParting(pair2)
        generateSkeletonSegment(p2,[x,y])
        
        if Segments[pSeg].previd==nSeg or Segments[nSeg].nextid==pSeg:
            x0 = x; y0 = y
            pair3=[pair2[1],pair1[0]]
            p3 = findEnd(pair3)
            generateSkeletonSegment(p3,[x,y])            
            continue
    
        tpSeg = Segments[pSeg].previd; tnSeg = Segments[nSeg].nextid
        if tpSeg==tnSeg:
            # Three or more sides and next segid = prior segid
            x0 = x; y0 = y
            x,y = intersectionBisectBisect(nSeg,tpSeg,tpSeg,pSeg)
            generateSkeletonSegment([x0,y0],[x,y])             
            pair1 = [nSeg,tpSeg]; pair2 = [tpSeg,pSeg]
            p1 = findEnd(pair1)   
            if p1 == None:
                ray = generateRay(Segments[nSeg],Segments[tpSeg])
                p1 = ray.p
                p1 = checkParting(pair1)
            generateSkeletonSegment(p1,[x,y])
        
            p2 = findEnd(pair2)
            if p2 != None:
                generateSkeletonSegment(p2,[x,y])
    return
#___________________________________________________________
title = ""

# Set up the bounds
# various Frillies
#bounds = (0, -600, 800, 0)

# Berg and star
bounds = (0, -1200, 1000, 0)

#CGal 2nd example
#bounds = (5, 45.5, 11, 48)

# Insert line Segments
Segments, numSegments = generateSegments()

# determine the corner bisectors and tracks
Corners, Tracks, Segments = generateTracks(Segments)
    
# Setup the graph
fig, ax = plt.subplots(figsize=(8, 8))
#___________________________________________________________

# Generate the motorcycle graph
# Determine intersection of reflex ray and closest corner
for track in Tracks:
    for i in range(numSegments):
        segment = Segments[i]
        cross, xi, yi, target =  \
         intersectionSegmentTrack(segment, track)
        if cross:
            track.bisect.p2[0]=xi
            track.bisect.p2[1]=yi            
            track.target = target
    
# Generate the intersection time pairs        
TIs = []
for track1 in Tracks:
    for track2 in Tracks:
        if track1.sid > track2.sid:
            cross,xi,yi = intersectionTrackTrack(track1,track2)            
            if cross:
                d1 = math.sqrt((xi-track1.bisect.p1[0])**2 + \
                           (yi-track1.bisect.p1[1])**2)
                d2 = math.sqrt((xi-track2.bisect.p1[0])**2 + \
                           (yi-track2.bisect.p1[1])**2)
                t1 = d1 * track1.speed
                t2 = d2 * track2.speed
                i = len(TIs)+1
                if t1<=t2:
                    TIs.append(TI(idsml=track1.tid, idlrg=track2.tid,\
                        tsml=t1, tlrg=t2, xi=xi, yi=yi))
                else:
                    TIs.append(TI(idsml=track2.tid, idlrg=track1.tid, \
                        tsml=t2, tlrg=t1, xi=xi, yi=yi))
                    
print("Number of edges:", numSegments)
L = len(Tracks)                    
print("Candidate pairs:", int(L*(L+1)/2))
print("Pairs used:",len(TIs))
print("")

# Sort the time pairs                    
def get_tsml(x):
    return x.tsml
            
TIs.sort(key=get_tsml)

# Determine death location
for TI in TIs:
    if(TI.tlrg < Tracks[TI.idlrg].t and (TI.tsml<Tracks[TI.idsml].t or \
            Tracks[TI.idsml].t==0) or Tracks[TI.idlrg].t==0):
        cross,xi,yi =  cross,xi,yi = \
            intersectionTrackTrack(Tracks[TI.idsml],Tracks[TI.idlrg])
        if cross:
            Tracks[TI.idlrg].t = TI.tlrg
            Tracks[TI.idlrg].clipid = TI.idsml
            Tracks[TI.idlrg].bisect.p2[0] = TI.xi
            Tracks[TI.idlrg].bisect.p2[1] = TI.yi
#___________________________________________________________

#"""
# Generate the skeleton

Skeletons = []
Ends      = []  # Unfinished coordinates from region removals
Parts     = []  # Data on pariing tracks
# Determine the regions to be eliminated
Regions = generateRegions()

for region in Regions:
    if Tracks[region.endtid].clipid != None:    # omit if end track is through
        found = findIsland(region, False)
        if found:
            found = findIsland(region, True)
            
# Find true target
findTrueTarget()
                
for region in Regions:
    span = (region.endid-region.startid)%numSegments
    i = Regions.index(region)
    if span>=2:
        #n,xi,yi = findFirst(region)
        n = findFirst(region)
        startSeg = spineMarch(n,region)
    if span==1:
        startSeg = spikeRemoval(region.startid)

"""        
print("Ends")
for i in range(len(Ends)):
    print(Ends[i])
"""

processConvexAreas()
#"""
  
#"""   
# plot the tracks
for track in Tracks:
    x1, y1 = track.bisect.p1
    x2, y2 = track.bisect.p2
    ax.plot([x1, x2], [y1, y2], color = "blue", linewidth=0.75 )
    ax.plot([x2],[y2], 'o', color ="blue", ms=1.75)
    
# Identify tracks
for track in Tracks:
    if track.tid!=None:
        x1, y1 = track.bisect.p1
        x2, y2 = track.bisect.p2
        ax.text(0.5*(x1+x2),0.5*(y1+y2), str(track.tid), \
             color="black", fontsize=7)
#"""


#___________________________________________________________
    
# Identify the corners and edges
for segment in Segments:
    x1, y1 = segment.p1
    x2, y2 = segment.p2
    if title=="Berg" and False:
        if segment.sid>=11 and segment.sid<=12 and segment.sid<numSegments:
            ax.plot([x1, x2], [y1, y2], color="red", linestyle="dashed")
        elif segment.sid>=7 and segment.sid<=10 and segment.sid<numSegments:
            ax.plot([x1, x2], [y1, y2], color="red", linestyle="dotted")
        elif segment.sid>=3 and segment.sid<=5 and segment.sid<numSegments:
            ax.plot([x1, x2], [y1, y2], color="red", linestyle="dashed")
        elif segment.sid>=13 and segment.sid<=14 and segment.sid<numSegments:            
            ax.plot([x1, x2], [y1, y2], color="red", linestyle="dotted")
        elif segment.sid==0:
            ax.plot([x1, x2], [y1, y2], color="red", linestyle="dotted")
        elif segment.sid==1 or segment.sid==2:
             ax.plot([x1, x2], [y1, y2], color="red", linestyle="dashdot")           
        elif segment.sid==6:
            ax.plot([x1, x2], [y1, y2], color="red")
    else:
        ax.plot([x1, x2], [y1, y2], color="red")        
        
    if segment.sid<numSegments:
        ax.text(0.5*(x1+x2),0.5*(y1+y2), str(segment.sid), \
         color="black", fontsize=7)
        ax.text(x1,y1, str(segment.sid), \
        color="red", fontsize=7)
    
ax.set_aspect('equal', 'box')
ax.set_xlim(bounds[0], bounds[2])
ax.set_ylim(bounds[1], bounds[3])
plt.title(title)
plt.xlabel("X")
plt.ylabel("Y")
plt.show()
     
