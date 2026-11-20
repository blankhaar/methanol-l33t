"""Streaming LAMDA checks, separate from the production serializer.

Exact reference hashes are checked by the release CLI. This reader also rejects
malformed/nonphysical serialization without loading millions of rate objects.
"""
import math
import re


def finite(token, *, positive=False):
    x = float(token)
    if not math.isfinite(x) or x < 0 or (positive and x == 0):
        raise ValueError('nonfinite or nonphysical number')
    return x


def identity(label):
    fields = label.split(':')
    sp, manifold = fields[:2]
    rest = fields[2:] if manifold.startswith('vt') else fields[3:]
    if sp not in ('A', 'E') or not re.fullmatch(r'J\d+', rest[0]) or not re.fullmatch(r'K[+-]?\d+', rest[1]):
        raise ValueError('invalid quantum label')
    j, k = int(rest[0][1:]), int(rest[1][1:])
    if abs(k) > j:
        raise ValueError('invalid J,K')
    p = None
    if sp == 'A':
        if k < 0 or len(rest) != 3 or rest[2] not in ('+', '-'):
            raise ValueError('invalid A parity')
        vt = int(manifold[2:]) if manifold.startswith('vt') else 0
        p = (1 if rest[2] == '+' else -1)*(-1)**(j+vt)
        if manifold == 'vnu10':
            p = -p
    elif len(rest) != 2:
        raise ValueError('invalid E label')
    return sp, j, p


def inspect_lamda(path):
    with open(path) as stream:
        def line():
            for raw in stream:
                if raw.strip() and not raw.lstrip().startswith(('!', '#')):
                    return raw.strip()
            raise ValueError('truncated LAMDA file')
        molecule = line()
        finite(line(), positive=True)
        nlev = int(line())
        if nlev < 2:
            raise ValueError('incomplete level system')
        levels = {}
        labels = set()
        previous = -1.0
        species = None
        for index in range(1, nlev+1):
            fields = line().split()
            if len(fields) != 4 or int(fields[0]) != index:
                raise ValueError('level index/width mismatch')
            energy, weight = finite(fields[1]), finite(fields[2], positive=True)
            label = fields[3]
            sp, j, p = identity(label)
            if energy <= previous or label in labels or weight != 2*j+1:
                raise ValueError('invalid level order, repeated identity or weight')
            if species is None:
                species = sp
            if sp != species:
                raise ValueError('A/E mixing')
            levels[index] = (energy, j, p)
            labels.add(label)
            previous = energy
        def endpoints(fields):
            u,l = int(fields[1]),int(fields[2])
            if u not in levels or l not in levels or u <= l:
                raise ValueError('missing or reversed endpoints')
            return u,l
        nrad = int(line())
        if nrad < 1:
            raise ValueError('missing radiative network')
        seen = set()
        adjacency = {i: set() for i in levels}
        for index in range(1, nrad+1):
            fields=line().split()
            if len(fields) != 6 or int(fields[0]) != index:
                raise ValueError('radiative index/width mismatch')
            u,l=endpoints(fields)
            ju,jl=levels[u][1],levels[l][1]
            if (u,l) in seen or abs(ju-jl)>1 or ju==jl==0:
                raise ValueError('duplicate or forbidden rank-one line')
            if species=='A' and levels[u][2]==levels[l][2]:
                raise ValueError('forbidden A parity line')
            finite(fields[3],positive=True);finite(fields[4],positive=True);finite(fields[5])
            seen.add((u,l));adjacency[u].add(l);adjacency[l].add(u)
        connected={1};frontier=[1]
        while frontier:
            i=frontier.pop()
            for neighbour in adjacency[i]-connected:
                connected.add(neighbour);frontier.append(neighbour)
        if len(connected)!=nlev:
            raise ValueError('disconnected radiative network')
        npart=int(line()); counts={};grids={}
        if npart!=3:
            raise ValueError('L33T requires three separate collider blocks')
        for _ in range(npart):
            collider=int(line().split()[0])
            if collider not in (2,3,6) or collider in counts:
                raise ValueError('wrong or duplicate collider')
            ncoll=int(line());nt=int(line())
            temps=[finite(x,positive=True) for x in line().split()]
            if len(temps)!=nt or any(a>=b for a,b in zip(temps,temps[1:])):
                raise ValueError('invalid temperature grid')
            if temps != [10,20,30,50,75,100,150,200,300,400]:
                raise ValueError('unexpected L33T grid')
            seen=set()
            for index in range(1,ncoll+1):
                fields=line().split()
                if len(fields)!=nt+3 or int(fields[0])!=index:
                    raise ValueError('collision index/width mismatch')
                pair=endpoints(fields)
                if pair in seen:
                    raise ValueError('duplicate collision pair')
                seen.add(pair)
                for token in fields[3:]:
                    finite(token)
            counts[collider]=ncoll;grids[collider]=temps
        if any(raw.strip() and not raw.lstrip().startswith(('!','#')) for raw in stream):
            raise ValueError('unexpected trailing records')
    return dict(levels=nlev,radiative=nrad,collision_rows=sum(counts.values()),
                collision_rows_by_partner=counts,temperature_grids_K=grids,
                radiative_connected=True)
