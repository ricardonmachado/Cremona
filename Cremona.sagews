︠f2551b5a-bbf6-4692-92df-cebf15339463︠
#log_matrix_list=i_M_n_d
︡38d5c0a5-fd25-4c02-a2e1-b17e7c7d96c7︡{"done":true}︡
︠7604ce79-7d44-4f21-be89-daff54ac3c66s︠
def natural_action(lista, permut):
    
    '''
    Def 2.13 para matrizes (apply_group_permut_action) (natural_action)
    Descrever melhor
    
    '''
    
    listatemp=[]
    
    for i in lista:
        listatemp.append([])
    
    for i in range(len(lista)):
        listatemp[i]=lista[permut[i]]

    return listatemp
︡c320a87e-ec8d-457a-80f9-fbeb0a6b11f9︡{"done":true}
︠55604b44-15a1-4d53-aa17-96e506444050s︠
def mat2list(mat):
    '''
    Transforma o objeto matriz no objeto lista, mantendo os dados
    
    '''
    newmat=[]
    nrow = mat.nrows()
    
    for i in range( nrow ):
        newmat.append( mat.rows()[i].list() )
    
    return newmat
︡ab322b54-389b-4f4c-9974-6e923a4c6dec︡{"done":true}
︠664f61a7-7714-4b02-8192-4d50ea4ad4e7s︠
def matrix_stabilizer(n,r, i_M_n_d):
    
    '''
    Pega cada matriz e gera a lista das permutações que não mudam as matrizes,
    no sentido de cremonas. (Estabilizador?)
    
    '''
    
    newlist = i_M_n_d[0]
    G=SymmetricGroup(n);
    Glist=G.list()
    GlistM=[]
    for g in Glist:
        GlistM.append(g.matrix())
    
    newlistG = []
    newlistG2 = []
    listaestabilizadora3=[]
    for nl in range(len(newlist)):
        listtemp=[]
        listtemp2=[]
        listaestabilizadora3.append([])
        for j in range(len(newlist[nl])):
            for k in range(len(mat2list(Matrix(newlist[nl][j])))):
                listtemp.append( Sequence( mat2list(Matrix(newlist[nl][j]))[k], immutable=True) ) #matriz como conjuntos
        listtemp = Set(listtemp)
        for g in range(len(Glist)):
            listtemp2.append([])
            for j in range(len(newlist[nl])):
                for k in range(len(mat2list(Matrix(newlist[nl][j])))):
                    listtemp2[g].append( Sequence(mat2list( (Glist[g].matrix()*(Matrix(newlist[nl][j]).transpose())).transpose() )[k], immutable=True) ) # lista do prod das matrizes como conj
            listtemp2[g] = Set(listtemp2[g])
        newlistG.append([])
        for g in range(len(Glist)):
            if listtemp2[g] == listtemp:
                newlistG[nl].append(Glist[g].matrix())
        
    for i in newlistG:
        if i != []:
            newlistG2.append(i)
    vetorbasico = vector(range(n))

    listaqqc=[]
    for i in range(len(newlistG2)):
        listaqqc.append([])
        listaqqc[i] = [newlistG2[i][k]*vetorbasico for k in range(len(newlistG2[i]))]

    return listaqqc
︡fbb219ef-5542-4822-b6d3-450a3c9df526︡{"done":true}
︠5e0f424d-eb2e-40d8-81bc-cb9f40e50b81s︠
def lemma_2_15(n, d, r, i_M_n_d, p):

    '''
    usamos o lemma 2.15 para eliminar casos onde
    a adição de uma nova coluna a uma determinada
    matriz log não gera uma órbita nova
    
    Para cada matrix newlist[i] da lista newlist e para cada j calcule a orbita de p2[i][j] via o estabilizador de newlist[i] via a acão de S_n e remova cada orbita da lista de candidatos a nova culuna de newlist[i]
    '''
    p2 = []
    newlist = i_M_n_d
    np = matrix_stabilizer(n,r, [newlist])

    for i in range(len(newlist)):
        
        p2.append(list(p)) # p2 é formado pelos candidatos a novas colunas
        #print 'p1', p2
        j = -1
        while j < len(p2[i])-1: # p2 é a imagem da ação de Sn sob f1. (ex. f1 =[1,1,1,0,0])
            j = j+1
            if list(p2[i][j]) not in newlist[i]:
                pjset = []
                for l in range(len(np[i])):
                    temp55 = list(p2[i][j])
                    temp = list(p2[i][j])
                    temp = natural_action(temp, np[i][l])
                    if temp != temp55 and temp not in pjset:  #Lema 2.15
                        pjset.append(temp) #  
                    
                    #(ini 1.3)#(p2[i] - pjset)###########################
                l2 = -1
                while l2 < len(p2[i])-1:
                    l2 = l2+1
                    temp22 = list(p2[i][l2])
                        
                    if temp22 in pjset:
                        p2[i].remove(temp22)
                            #(fim 1.3)############################
    #print 'p2', p2                        
    return p2

︡22f04216-cee6-4f05-a9fb-d490d6d03cd5︡{"done":true}
︠7b28df7e-cb84-477c-824d-7f5f025f5751s︠
#(esta rotina acrescenta mais uma coluna nos representantes das orbitas)######
def add_col_lemma_2_14(d, p2, newlist):
    '''
    We also use Proposition 2.6 in lines 21 to 30
    
    '''
    
    newlist3 = []
    newlist2 = []
    for i in range(len(newlist)):
        for j in range(len(p2[i])):

            if p2[i][j] in newlist[i]:
                newlist2 = []
                
            else:
                tempx = list(newlist[i]) # faz copia 
                tempx.append(p2[i][j])
                newlist2 = tempx

            if len(newlist2)>0:
                m = Matrix(newlist2)
                #print m.echelon_form()
                if m.echelon_form()[-1] != 0:

                    if m.nrows() != m.ncols():
                        newlist3.append(newlist2)
                    else:
                        if abs(m.det()) == d:
                            newlist3.append(newlist2)
                            
    return newlist3
︡0732c970-120b-47cc-853f-201b44feae87︡{"done":true}
︠2c4016eb-7ba5-406c-9a9f-4dc3cf5df947s︠
def i_M_n_d_by_incidence_degree(n,lista): #Lemma_2_18
    
    '''
    Separa a i_M_n_d em classes, na qual, cada classe é uma lista de graus de incidência
    '''
    
    grupoi=vector([0]*n)
    par=[]
    grupo=[]
    
    for i in range(len(lista)):
        grupo.append(grupoi)
        for j in range(len(lista[0])):
            grupo[i] = grupo[i]+vector(lista[i][j])
        
        a = list(grupo[i])
        a.sort()
        grupo[i]=a
        par.append([grupo[i], lista[i]])
    
    degree_list=[]
    for i in range(len(grupo)):
        a = list(grupo[i])
        a.sort()
        degree_list.append(a)
    
    k = degree_list
    import itertools
    k.sort()
    degree_list = list(k for k,_ in itertools.groupby(k))
    
    degree_bag=[]
    new_bag=[]
    for i in range(len(degree_list)):
        new_bag.append([])
        degree_bag.append(degree_list[i])
        for j in range(len(par)):
            if par[j][0] == degree_list[i]:
                new_bag[i].append(par[j][1])

    return new_bag
︡9ac95f05-5f52-42b2-b3cb-4d5271075ac3︡{"done":true}
︠c43d7507-0aea-4599-b2d7-74969e51a69as︠
#Aplica o Lema 2.14 e 2.15, 4
def build_next_i_M_n_d(n,d,r, i_M_n_d):
    '''
    
    Aplica os lemas 2.14 e 2.15
    
    '''
    
    temp_i_M_n_d = []
    
    temp1 = [1]*d
    temp0 = [0]*(n-d)
    primary_log_vector = temp1+temp0 # primary_log_vector=f1
    p = Permutations(primary_log_vector).list()
    
    for i in range(len(p)):
        p[i] = list(p[i])
    
    if r == 1:
        temp_i_M_n_d = [[primary_log_vector]]
        return [temp_i_M_n_d]
    
    else: #if r >= 2:
        newlist = i_M_n_d[0]
        
        p2 = lemma_2_15(n, d, r, newlist, p)
        
        temp_i_M_n_d = add_col_lemma_2_14(d, p2, newlist)
        
        temp_i_M_n_d = i_M_n_d_by_incidence_degree(n,temp_i_M_n_d)
        
        
        return temp_i_M_n_d
︡de3acf84-ff1d-4054-8591-2c4edfcf3485︡{"done":true}
︠0a48ff43-bb7a-4843-ad6d-3fd54780c72cs︠
def orbits(n,d, fiset):
    
    '''
    Gera a orbita o primeiro elemento da lista fiset do input
    '''
    
    orb1=[]
    orb2=[]
    fiset1 = fiset[0]
    fiset2 = fiset[0]    
    
    permut1 = Permutations(range(n))
    m = len(fiset1)
    permut2 = Permutations(range(m))
    
    for i in range(len(permut2)):
        orb2.append(natural_action(fiset2, permut2[i] ))
    
    for j in range(len(orb2)):
        orb2[j] =  Matrix(orb2[j]).transpose()
        orb2[j] = mat2list(orb2[j])
        
        for i in range(len(permut1)):
            orb1.append( natural_action(orb2[j], permut1[i] ) )
            #orb1.append( Sequence( apllypermutation2(orb2[j], permut1[i] ) , immutable=True) )###############
    
    k = orb1                                          #
    import itertools                                  #
    k.sort()                                          # sort
    orb1 = list(k for k,_ in itertools.groupby(k))    #
    
    return orb1
︡f452b69b-bd18-4d94-a9bd-39ff60b73e0f︡{"done":true}
︠eac84de4-1ba9-432d-82af-c14811723ec9s︠
# talvez use o Lema 2.15
def refresh_temp_i_M_n_d(n,d, fiset):
    
    '''
    Dada a lista de matrizes log, essa rotina elimina dessa lista, 
    a orbita do seu primeiro elemento
    
    '''
    #print fiset
    orb = orbits(n,d, fiset)
    x=len(fiset)
    l=range(x)
    
    newset1 = [mat2list( Matrix(fiset[i]).transpose() ) for i in l]
    
    #a = time()
    newset2 = [newset1[i] for i in l if (newset1[i] not in orb)] #newset1 - orb
    #print time()-a
    
    l2 = range(len(newset2))
    teste = [mat2list( Matrix(newset2[i]).transpose() ) for i in l2]
    
    return teste
︡b5a7eaa0-163b-4b7f-a4f8-3e9e320b881b︡{"done":true}
︠b8fae6d2-ad33-455c-8349-0168e2c65d2ds︠
def final_list(n,d):
    
    '''
    Rotina do fluxograma
    
    '''
    
    temp_i_M_n_d=[] 
    i_M_n_d = []
    print '(n_col, n_it, n_class, n_elts_class)'
    
    for i in range(1,n+1):
        set = []
        temp_i_M_n_d = build_next_i_M_n_d(n, d, i, i_M_n_d)
        i_M_n_d=[]
        k=0
        for j in range(len(temp_i_M_n_d)):
            #print 'teste'
            temp = temp_i_M_n_d[j]
            while len(temp)>0:
                k=k+1
                print '(   %d' %i,',  %d' %k,',     %d' %len(temp_i_M_n_d),',        %d' %len(temp_i_M_n_d[j]), '  )'
                
                set.append(temp[0])
                temp = refresh_temp_i_M_n_d(n,d, temp)
        
        i_M_n_d.append(set)
    #log_matrix_list=1    
    #listafinalc = []
    #m = set
    #set=0    
    
    #for i in range(len(m)):
    #    if Matrix(m[i]).det() == d or Matrix(m[i]).det() == -d:
    #        listafinalc.append(m[i])
    
    #return [len(listafinalc),listafinalc]
    return [len(set), set]
    
︡849ceea4-f5b0-4caf-8698-2d63abf5aa9d︡{"done":true}
︠d133a7f5-b4ad-4080-b254-da11f8d57a43s︠
final_list(5,3)
︡d4aa45a9-718c-4d52-be8f-9a597a5f051b︡{"stdout":"(n_col, n_it, n_class, n_elts_class)\n(   1 ,  1 ,     1 ,        1   )\n"}︡{"stdout":"(   2"}︡{"stdout":" ,  1 ,     2 ,        1   )\n(   2 ,  2 ,     2 ,        1   )\n(   3"}︡{"stdout":" ,  1 ,     4 ,        1   )\n(   3"}︡{"stdout":" ,  2 ,     4 ,        1   )\n(   3"}︡{"stdout":" ,  3 ,     4 ,        3   )\n(   3"}︡{"stdout":" ,  4 ,     4 ,        2   )\n(   4"}︡{"stdout":" ,  1 ,     4 ,        1   )\n(   4"}︡{"stdout":" ,  2 ,     4 ,        3   )\n(   4"}︡{"stdout":" ,  3 ,     4 ,        3   )\n(   4"}︡{"stdout":" ,  4 ,     4 ,        4   )\n(   4"}︡{"stdout":" ,  5 ,     4 ,        4   )\n(   5"}︡{"stdout":" ,  1 ,     4 ,        3   )\n(   5"}︡{"stdout":" ,  2 ,     4 ,        3   )\n(   5"}︡{"stdout":" ,  3 ,     4 ,        4   )\n(   5"}︡{"stdout":" ,  4 ,     4 ,        1   )\n[4, [[[1, 1, 1, 0, 0], [1, 0, 1, 1, 0], [1, 1, 0, 1, 0], [0, 1, 1, 1, 0], [1, 0, 0, 1, 1]], [[1, 1, 1, 0, 0], [1, 0, 1, 1, 0], [1, 1, 0, 1, 0], [1, 0, 0, 1, 1], [0, 1, 0, 1, 1]], [[1, 1, 1, 0, 0], [1, 0, 1, 1, 0], [1, 1, 0, 1, 0], [1, 0, 0, 1, 1], [0, 1, 1, 0, 1]], [[1, 1, 1, 0, 0], [1, 0, 1, 1, 0], [0, 1, 1, 0, 1], [1, 0, 0, 1, 1], [0, 1, 0, 1, 1]]]]"}︡{"stdout":"\n"}︡{"done":true}
︠063a26ef-48e4-4dd9-88bc-e53d4fc9e988︠

︡f76fb7b5-a866-4f23-8444-f18f54192738︡
︠bf3f7c2f-e7f9-460c-87a3-389476238506s︠
temp= final_list(5,3)
︡00ad0f4b-c59b-4af6-96c5-2ea8395587b0︡{"stdout":"(n_col, n_it, n_class, n_elts_class)\n(   1 ,  1 ,     1 ,        1   )\n(   2"}︡{"stdout":" ,  1 ,     2 ,        1   )\n(   2 ,  2 ,     2 ,        1   )\n(   3"}︡{"stdout":" ,  1 ,     4 ,        1   )\n(   3"}︡{"stdout":" ,  2 ,     4 ,        1   )\n(   3"}︡{"stdout":" ,  3 ,     4 ,        3   )\n(   3"}︡{"stdout":" ,  4 ,     4 ,        2   )\n(   4"}︡{"stdout":" ,  1 ,     4 ,        1   )\n(   4"}︡{"stdout":" ,  2 ,     4 ,        3   )\n(   4"}︡{"stdout":" ,  3 ,     4 ,        3   )\n(   4"}︡{"stdout":" ,  4 ,     4 ,        4   )\n(   4"}︡{"stdout":" ,  5 ,     4 ,        4   )\n(   5"}︡{"stdout":" ,  1 ,     4 ,        3   )\n(   5"}︡{"stdout":" ,  2 ,     4 ,        3   )\n(   5"}︡{"stdout":" ,  3 ,     4 ,        4   )\n(   5"}︡{"stdout":" ,  4 ,     4 ,        1   )\n"}︡{"done":true}︡
︠71e9b3b7-4e8e-41c0-b12b-17a2e1ed30d7s︠
def list_2_mono(lista):
        
    n = len(lista)
    
    R, x = PolynomialRing(ZZ, 'x', n+1).objgens()
    x = x[1:n+1]
    temp = range(1,n+1)
    
    temp2=[]
    for j in range(n):
        temp2.append([])
        for i in range(n):
            temp2[j].append(x[i]*lista[j][i])

            
    mon_list=[]
    for j in temp2:
        mon = 1
        for i in j:
            if i != 0:
                mon = mon*i
        mon_list.append(mon)
        
    
    return mon_list
︡d0a23ad9-300a-4d43-bb44-41e323c32843︡{"done":true}︡
︠f9713f20-76c2-48b0-be97-7c0f7ec23b41s︠
def list_2_mono_final(lista):
    
    n = len(lista)
    
    lista_final=[]
    for i in lista:
        lista_final.append(list_2_mono(i))
    
    return lista_final
︡a84886b3-ff0b-4a8f-9d4d-828dfd140df4︡{"done":true}︡
︠ac9526eb-3671-4f72-a8e7-c057b214e14es︠
temp[1]
︡aba48bc7-6b2e-4658-b844-22aeaa2b87c2︡{"stderr":"Error in lines 1-1\nTraceback (most recent call last):\n  File \"/projects/sage/sage-7.5/local/lib/python2.7/site-packages/smc_sagews/sage_server.py\", line 982, in execute\n    exec compile(block+'\\n', '', 'single') in namespace, locals\n  File \"\", line 1, in <module>\nNameError: name 'temp' is not defined\n"}︡{"done":true}︡
︠8cb59994-aec8-4897-8e3b-1111ae6f241as︠
list_2_mono_final(temp[1])
︡941a934b-71df-4bb9-a4cc-becb22d32427︡
︠1bd9ac8d-c9b8-4464-8665-dbeb958604d0︠
final_list(6,2)
︡5367a78d-3fef-4179-8d69-8e5da3f2191c︡{"stdout":"(n_col, n_it, n_class, n_elts_class)\n(   1 ,  1 ,     1 ,        1   )\n(   2"}︡{"stdout":" ,  1 ,     2 ,        2   )\n(   2"}︡{"stdout":" ,  2 ,     2 ,        1   )\n(   3"}︡{"stdout":" ,  1 ,     5 ,        1   )\n(   3"}︡{"stdout":" ,  2 ,     5 ,        1   )\n(   3"}︡{"stdout":" ,  3 ,     5 ,        2   )\n(   3"}︡{"stdout":" ,  4 ,     5 ,        2   )\n(   3"}︡{"stdout":" ,  5 ,     5 ,        1   )\n(   4"}︡{"stdout":" ,  1 ,     6 ,        3   )\n(   4"}︡{"stdout":" ,  2 ,     6 ,        1   )\n(   4"}︡{"stdout":" ,  3 ,     6 ,        3   )\n(   4"}︡{"stdout":" ,  4 ,     6 ,        4   )\n(   4"}︡{"stdout":" ,  5 ,     6 ,        4   )\n(   4"}︡{"stdout":" ,  6 ,     6 ,        2   )\n(   4"}︡{"stdout":" ,  7 ,     6 ,        5   )\n(   4"}︡{"stdout":" ,  8 ,     6 ,        5   )\n(   5"}︡{"stdout":" ,  1 ,     9 ,        3   )\n(   5"}︡{"stdout":" ,  2 ,     9 ,        3   )\n(   5"}︡{"stdout":" ,  3 ,     9 ,        4   )\n(   5"}︡{"stdout":" ,  4 ,     9 ,        1   )\n(   5"}︡{"stdout":" ,  5 ,     9 ,        1   )\n(   5"}︡{"stdout":" ,  6 ,     9 ,        4   )\n(   5"}︡{"stdout":" ,  7 ,     9 ,        2   )\n(   5"}︡{"stdout":" ,  8 ,     9 ,        11   )\n(   5"}︡{"stdout":" ,  9 ,     9 ,        11   )\n(   5"}︡{"stdout":" ,  10 ,     9 ,        11   )\n(   5"}︡{"stdout":" ,  11 ,     9 ,        5   )\n(   5"}︡{"stdout":" ,  12 ,     9 ,        5   )\n(   6"}︡{"stdout":" ,  1 ,     6 ,        4   )\n(   6"}︡{"stdout":" "}︡{"stdout":",  2 ,     6 ,        5   )\n(   6"}︡{"stdout":" "}︡{"stdout":",  3 ,     6 ,        2   )\n(   6"}︡{"stdout":" "}︡{"stdout":",  4 ,     6 ,        5   )\n(   6"}︡{"stdout":" "}{"stdout":" ,  6 ,     9 ,        4   )\n(   5"}︡{"stdout":" ,  7 ,     9 ,        2   )\n(   5"}︡{"stdout":" ,  8 ,     9 ,        11   )\n(   5"}︡{"stdout":" ,  9 ,     9 ,        11   )\n(   5"}︡{"stdout":" ,  10 ,     9 ,        11   )\n(   5"}︡{"stdout":" ,  11 ,     9 ,        5   )\n(   5"}︡{"stdout":" ,  12 ,     9 ,        5   )\n(   6"}︡{"stdout":" ,  1 ,     6 ,        4   )\n(   6"}︡{"stdout":" "}︡{"stdout":",  2 ,     6 ,        5   )\n(   6"}︡{"stdout":" "}︡{"stdout":",  3 ,     6 ,        2   )\n(   6"}︡{"stdout":" "}︡{"stdout":",  4 ,     6 ,        5   )\n(   6"}︡{"stdout":" "}︡
︠cbb56f67-c51b-4967-804d-66c4521460c4︠









