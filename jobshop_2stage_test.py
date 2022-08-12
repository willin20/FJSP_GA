import matplotlib.pyplot as plt
import numpy as np
import pprint
import random as r
import copy
import time

#交叉
def pox(p1_os,p2_os,p1_ma,p2_ma,N_jobs,N_machine,sum_opr):       #(親個体1,2(操作)、親個体1,2(機械),仕事数,機械数,操作数)
    c1_os = [N_jobs]*sum_opr        #子個体１(操作)
    c2_os = [N_jobs]*sum_opr         #子個体２(操作)
    c1_ma = [N_machine]*sum_opr      #子個体1(機械)
    c2_ma = [N_machine]*sum_opr     #子個体2(機械)
    sub_1 = ([])                    #小集団1
    sub_2 = ([])                    #小集団２
    for i in range(N_jobs):         #仕事を２つの集合に分ける
        if r.random()<0.5:
            sub_1.append(i)
        else:
            sub_2.append(i)
    for i in range(sum_opr):             #親１にある集合１の仕事は子1の同じ場所へ
        if p1_os[i] in sub_1:
            c1_os[i] = p1_os[i]
            c1_ma[i] = p1_ma[i]
    for i in range(sum_opr):             #親2にある集合1の仕事は子2の同じ場所へ
        if p2_os[i] in sub_1:
            c2_os[i] = p2_os[i]
            c2_ma[i] = p2_ma[i]
    for i in range(sum_opr):             #親１にある集合2の仕事は子2へ順番を変えずに引き継ぐ
        if p1_os[i] in sub_2:
            if N_jobs in c2_os:
                a = c2_os.index(N_jobs)
                c2_os[a] = p1_os[i]
                c2_ma[a] = p1_ma[i]
    for i in range(sum_opr):             #親2にある集合2の仕事は子1へ順番を変えずに引き継ぐ
        if p2_os[i] in sub_2:
            if N_jobs in c1_os:
                a = c1_os.index(N_jobs)
                c1_os[a] = p2_os[i]
                c1_ma[a] = p2_ma[i]
    return c1_os,c2_os,c1_ma,c2_ma

#突然変異(操作)
def swap_mutation_os(gene_os,gene_ma,N):          #2点間の入れ替え #(遺伝子、遺伝子数)
    m_i = r.randint(0,N-1)
    m_j = r.randint(0,N-1)
    gene_os[m_i],gene_os[m_j] = gene_os[m_j],gene_os[m_i]
    gene_ma[m_i],gene_ma[m_j] = gene_ma[m_j],gene_ma[m_i]
    return gene_os,gene_ma

#突然変異(機械)
def mutation_ma(gene_ma,gene_os,N_machine,num_opr,opr_list,mat_ope):
    opr_remain_list = copy.copy(num_opr)
    machine_list = [i for i in range(N_machine)]
    for i in range(sum_opr):
        p_sel_machine = []                                                          #操作する機械の選択確率
        for j in range(N_machine):
            if mat_ope[opr_list[gene_os[i]]+num_opr[gene_os[i]]-opr_remain_list[gene_os[i]]][j] != 10**5:
                p_sel_machine.append(mat_ope[opr_list[gene_os[i]]+num_opr[gene_os[i]]-opr_remain_list[gene_os[i]]][j])
            else:
                p_sel_machine.append(0)
        s = sum(p_sel_machine)
        if sum([i != 0 for i in p_sel_machine]) == 1:
            gene_ma[i] = p_sel_machine.index(max(p_sel_machine))
        else:
            for j in range(N_machine):
                if mat_ope[opr_list[gene_os[i]]+num_opr[gene_os[i]]-opr_remain_list[gene_os[i]]][j] != 10**5:
                    p_sel_machine[j] = s-mat_ope[opr_list[gene_os[i]]+num_opr[gene_os[i]]-opr_remain_list[gene_os[i]]][j]
            sel_machine = r.choices(machine_list,weights=p_sel_machine)
            gene_ma[i] = sel_machine[0]
        opr_remain_list[gene_os[i]] = opr_remain_list[gene_os[i]]-1
    return gene_ma
    
#初期化    
def initialization(num_opr,N_machine,opr_list,mat_ope,N_job):     #(操作数、機械数、操作数の区切り)
    opr_remain_list = copy.copy(num_opr)
    machine_list = [i for i in range(N_machine)]
    sol_MA = []
    sol_OS = []
    while opr_remain_list != [0]*N_job:
        job_max_list = [i for i, v in enumerate(opr_remain_list) if v == max(opr_remain_list)]      #残り操作数が最大の仕事
        job_sel = r.choice(job_max_list)                                            #選んだ仕事
        sol_OS.append(job_sel)
        p_sel_machine = []                                                          #操作する機械の選択確率
        for i in range(N_machine):
            if mat_ope[opr_list[job_sel]+num_opr[job_sel]-opr_remain_list[job_sel]][i] != 10**5:
                # p_sel_machine.append(mat_ope[opr_list[job_sel]+num_opr[job_sel]-opr_remain_list[job_sel]][i]/sum(mat_ope[opr_list[job_sel]+num_opr[job_sel]-opr_remain_list[job_sel]]))
                p_sel_machine.append(mat_ope[opr_list[job_sel]+num_opr[job_sel]-opr_remain_list[job_sel]][i])
            else:
                p_sel_machine.append(0)
        s = sum(p_sel_machine)
        # print(p_sel_machine)
        if sum([i != 0 for i in p_sel_machine]) == 1:
            sol_MA.append(p_sel_machine.index(max(p_sel_machine)))
        else:    
            for i in range(N_machine):
                if mat_ope[opr_list[job_sel]+num_opr[job_sel]-opr_remain_list[job_sel]][i] != 10**5:
                    # p_sel_machine.append(mat_ope[opr_list[job_sel]+num_opr[job_sel]-opr_remain_list[job_sel]][i]/sum(mat_ope[opr_list[job_sel]+num_opr[job_sel]-opr_remain_list[job_sel]]))
                    p_sel_machine[i] = s-mat_ope[opr_list[job_sel]+num_opr[job_sel]-opr_remain_list[job_sel]][i]
            sel_machine = r.choices(machine_list,weights=p_sel_machine)
            sol_MA.append(sel_machine[0])
        opr_remain_list[job_sel] = opr_remain_list[job_sel]-1
    return sol_OS,sol_MA

#選択
def selection(gene_os,gene_ma,evaluation,N_gene):
    list = [i for i in range(N_gene)]
    new_gene_os = []
    new_gene_ma = []
    new_gene_num = r.choices(list,k=N_gene,weights=evaluation)          #ルーレット選択
    for i in range(N_gene):                                             #選んだ個体を代入
        new_gene_os.append(gene_os[new_gene_num[i]])
        new_gene_ma.append(gene_ma[new_gene_num[i]])
    return new_gene_os,new_gene_ma

#エリート保存
def elitizm(gene_os,gene_ma,makespan,N_gene):
    list_elite = [10**5]*int(N_gene/10)             #上位10%の個体のメイクスパンを格納するリスト
    list_elite_index = [0]*int(N_gene/10)           #上位10%の個体のインデックスを格納するリスト
    gene_elite_os = []                              #上位10%の個体の機械割当てを格納するリスト
    gene_elite_ma = []                              #上位10%の個体の操作順序を格納するリスト
    for i in range(N_gene):
        if max(list_elite) > makespan[i]:
            list_elite_index[list_elite.index(max(list_elite))] = i                 #先にインデックスを変更、もしくは変更前のインデックスを保存
            list_elite[list_elite.index(max(list_elite))] = makespan[i]             #先に値を変更するとインデックスが変わるから
    for i in range(len(list_elite_index)):          #上位個体の機械割当て、操作順序のリストを代入
        gene_elite_os.append(gene_os[list_elite_index[i]])
        gene_elite_ma.append(gene_ma[list_elite_index[i]])
    return gene_elite_os,gene_elite_ma

#各個体におけるメイクスパンの評価
def makespan(opr_MA,opr_OS,num_machine,num_operation,benchmark,num_job,opr_list):   #(機械割当て遺伝子列、操作順序遺伝子列、機械数、操作数、工程表、仕事数、操作数の区切り)
    machine_available_time = [0]*num_machine                #各機械の利用可能時間
    job_available_time = [0]*num_job                        #各仕事の割当て可能時間
    operation_graph = [[0]*3 for i in range(num_operation)]   #操作開始時間、操作終了時間、割当て機械の順
    operaton_cnt = [0]*num_job                               #各仕事の完了した操作数カウンタ
    for i in range(num_operation):
        if machine_available_time[opr_MA[i]]>job_available_time[opr_OS[i]]:     #機械の操作可能開始時間>仕事の割当て開始可能時間
            operation_graph[i][0] = machine_available_time[opr_MA[i]]
            operation_graph[i][1] = machine_available_time[opr_MA[i]]+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph[i][2] = opr_MA[i]
        else:                                                                   #機械の操作可能開始時間<仕事の割当て開始可能時間
            operation_graph[i][0] = job_available_time[opr_OS[i]]
            operation_graph[i][1] = job_available_time[opr_OS[i]]+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph[i][2] = opr_MA[i]
        machine_available_time[opr_MA[i]] = operation_graph[i][1]
        job_available_time[opr_OS[i]] = operation_graph[i][1]
        operaton_cnt[opr_OS[i]] += 1
    return max(job_available_time)

#ガントチャート
def ganttchart(opr_MA,opr_OS,num_machine,num_operation,benchmark,num_job,opr_list):   #(機械割当て遺伝子列、操作順序遺伝子列、機械数、操作数、工程表、仕事数、操作数の区切り)
    machine_available_time = [0]*num_machine                #各機械の利用可能時間
    job_available_time = [0]*num_job                        #各仕事の割当て可能時間
    operation_graph = [[0]*3 for i in range(num_operation)]   #操作開始時間、操作終了時間、割当て機械の順
    operaton_cnt = [0]*num_job                               #各仕事の完了した操作数カウンタ
    cmap = plt.get_cmap("tab10")            #仕事ごとに色を決定
    for i in range(num_operation):
        if machine_available_time[opr_MA[i]]>job_available_time[opr_OS[i]]:     #機械の操作可能開始時間>仕事の割当て開始可能時間
            operation_graph[i][0] = machine_available_time[opr_MA[i]]
            operation_graph[i][1] = machine_available_time[opr_MA[i]]+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph[i][2] = opr_MA[i]
        else:                                                                   #機械の操作可能開始時間<仕事の割当て開始可能時間
            operation_graph[i][0] = job_available_time[opr_OS[i]]
            operation_graph[i][1] = job_available_time[opr_OS[i]]+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph[i][2] = opr_MA[i]
        machine_available_time[opr_MA[i]] = operation_graph[i][1]     #機械の操作可能時間の更新
        job_available_time[opr_OS[i]] = operation_graph[i][1]         #仕事の割当て可能時間の更新
        plt.barh(y=operation_graph[i][2],width=benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]],left=operation_graph[i][0],color=cmap(opr_OS[i]),ec='k')    #プロット
        plt.text((operation_graph[i][0]+operation_graph[i][1])/2,operation_graph[i][2],opr_OS[i],ha="center")
        operaton_cnt[opr_OS[i]] += 1        #各仕事の完了した操作数を1増やす
    plt.xlabel("Makespan")
    plt.ylabel("machine")
    plt.show()          #図の提示

#故障する機械の選択と選択された機械の負荷を出力
def select_break_machine(population_OS,population_MA,num_machine,num_operation,num_job,benchmark,opr_list):
    list_machine = [i for i in range(num_machine)]
    list_machine_workroad = [0]*num_machine     #機械にかかる負荷
    operation_cnt = [0]*num_job
    for i in range(num_operation):              #機械にかかる負荷を計算
        list_machine_workroad[population_MA[i]] += benchmark[opr_list[population_OS[i]]+operation_cnt[population_OS[i]]][population_MA[i]]
        operation_cnt[population_OS[i]] += 1
    break_machine = r.choices(list_machine,weights=list_machine_workroad)      #故障する機械の選択
    MBT = list_machine_workroad[break_machine[0]]
    return break_machine[0],MBT

#各個体におけるメイクスパンの評価(機械の故障)
def makespan_machine_break(opr_MA,opr_OS,num_machine,num_operation,benchmark,num_job,opr_list,tau_m,tau_s,tau_d):   #(機械割当て遺伝子列、操作順序遺伝子列、機械数、操作数、工程表、仕事数、操作数の区切り)
    machine_available_time = [0]*num_machine                #各機械の利用可能時間
    job_available_time = [0]*num_job                        #各仕事の割当て可能時間
    operation_graph = [[0]*3 for i in range(num_operation)]   #操作開始時間、操作終了時間、割当て機械の順
    operaton_cnt = [0]*num_job                               #各仕事の完了した操作数カウンタ
    for i in range(num_operation):
        if machine_available_time[opr_MA[i]]>job_available_time[opr_OS[i]]:     #機械の操作可能開始時間>仕事の割当て開始可能時間
            operation_graph[i][0] = machine_available_time[opr_MA[i]]
            operation_graph[i][1] = machine_available_time[opr_MA[i]]+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph[i][2] = opr_MA[i]
        else:                                                                   #機械の操作可能開始時間<仕事の割当て開始可能時間
            operation_graph[i][0] = job_available_time[opr_OS[i]]
            operation_graph[i][1] = job_available_time[opr_OS[i]]+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph[i][2] = opr_MA[i]
        if ((operation_graph[i][0] <= tau_s <= operation_graph[i][1]) or (operation_graph[i][0] <= tau_s+tau_d <= operation_graph[i][1])) and operation_graph[i][2] == tau_m:
            operation_graph[i][0] = tau_s+tau_d
            operation_graph[i][1] = tau_s+tau_d+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph[i][2] = opr_MA[i]
        machine_available_time[opr_MA[i]] = operation_graph[i][1]
        job_available_time[opr_OS[i]] = operation_graph[i][1]
        operaton_cnt[opr_OS[i]] += 1
    return max(job_available_time)

#故障後のメイクスパンと故障前後の操作完了時間のずれの合計を計算
def calculation_MS_and_STB(opr_MA,opr_OS,num_machine,num_operation,benchmark,num_job,opr_list,a,b):   #(機械割当て遺伝子列、操作順序遺伝子列、機械数、操作数、工程表、仕事数、操作数の区切り)
    #故障前のメイクスパンを計算
    machine_available_time = [0]*num_machine                #各機械の利用可能時間
    job_available_time = [0]*num_job                        #各仕事の割当て可能時間
    operation_graph_p = [[0]*3 for i in range(num_operation)]   #操作開始時間、操作終了時間、割当て機械の順
    operaton_cnt = [0]*num_job                               #各仕事の完了した操作数カウンタ
    for i in range(num_operation):      #各操作のガントチャートを作成
        if machine_available_time[opr_MA[i]]>job_available_time[opr_OS[i]]:     #機械の操作可能開始時間>仕事の割当て開始可能時間
            operation_graph_p[i][0] = machine_available_time[opr_MA[i]]
            operation_graph_p[i][1] = machine_available_time[opr_MA[i]]+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph_p[i][2] = opr_MA[i]
        else:                                                                   #機械の操作可能開始時間<仕事の割当て開始可能時間
            operation_graph_p[i][0] = job_available_time[opr_OS[i]]
            operation_graph_p[i][1] = job_available_time[opr_OS[i]]+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph_p[i][2] = opr_MA[i]
        machine_available_time[opr_MA[i]] = operation_graph_p[i][1]
        job_available_time[opr_OS[i]] = operation_graph_p[i][1]
        operaton_cnt[opr_OS[i]] += 1
    #故障後のメイクスパンを計算
    machine_available_time = [0]*num_machine                #各機械の利用可能時間
    job_available_time = [0]*num_job                        #各仕事の割当て可能時間
    operation_graph_q = [[0]*3 for i in range(num_operation)]   #操作開始時間、操作終了時間、割当て機械の順
    operaton_cnt = [0]*num_job                               #各仕事の完了した操作数カウンタ
    tau_m, MBT = select_break_machine(opr_OS,opr_MA,num_machine,num_operation,num_job,benchmark,opr_list)
    tau_s = a*MBT
    tau_d = b*MBT
    for i in range(num_operation):      #各操作のガントチャートを作成
        if machine_available_time[opr_MA[i]]>job_available_time[opr_OS[i]]:     #機械の操作可能開始時間>仕事の割当て開始可能時間
            operation_graph_q[i][0] = machine_available_time[opr_MA[i]]
            operation_graph_q[i][1] = machine_available_time[opr_MA[i]]+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph_q[i][2] = opr_MA[i]
        else:                                                                   #機械の操作可能開始時間<仕事の割当て開始可能時間
            operation_graph_q[i][0] = job_available_time[opr_OS[i]]
            operation_graph_q[i][1] = job_available_time[opr_OS[i]]+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph_q[i][2] = opr_MA[i]
        #故障期間と操作実行時間が被っている場合、機械復旧後に再割り当て
        if ((operation_graph_q[i][0] <= tau_s <= operation_graph_q[i][1]) or (operation_graph_q[i][0] <= tau_s+tau_d <= operation_graph_q[i][1])) and operation_graph_q[i][2] == tau_m:
            operation_graph_q[i][0] = tau_s+tau_d
            operation_graph_q[i][1] = tau_s+tau_d+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph_q[i][2] = opr_MA[i]
        machine_available_time[opr_MA[i]] = operation_graph_q[i][1]
        job_available_time[opr_OS[i]] = operation_graph_q[i][1]
        operaton_cnt[opr_OS[i]] += 1
    #故障前後の各操作の完了時間のずれを計算
    M = 0
    for i in range(num_operation):
        M += abs(operation_graph_p[i][1]-operation_graph_q[i][1])
    return max(job_available_time),M

#ガントチャート(機械の故障あり)
def ganttchart_machine_break(opr_MA,opr_OS,num_machine,num_operation,benchmark,num_job,opr_list,tau_m,tau_s,tau_d):   #(機械割当て遺伝子列、操作順序遺伝子列、機械数、操作数、工程表、仕事数、操作数の区切り)
    machine_available_time = [0]*num_machine                #各機械の利用可能時間
    job_available_time = [0]*num_job                        #各仕事の割当て可能時間
    operation_graph = [[0]*3 for i in range(num_operation)]   #操作開始時間、操作終了時間、割当て機械の順
    operaton_cnt = [0]*num_job                               #各仕事の完了した操作数カウンタ
    cmap = plt.get_cmap("tab10")            #仕事ごとに色を決定
    for i in range(num_operation):
        if machine_available_time[opr_MA[i]]>job_available_time[opr_OS[i]]:     #機械の操作可能開始時間>仕事の割当て開始可能時間
            operation_graph[i][0] = machine_available_time[opr_MA[i]]
            operation_graph[i][1] = machine_available_time[opr_MA[i]]+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph[i][2] = opr_MA[i]
        else:                                                                   #機械の操作可能開始時間<仕事の割当て開始可能時間
            operation_graph[i][0] = job_available_time[opr_OS[i]]
            operation_graph[i][1] = job_available_time[opr_OS[i]]+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph[i][2] = opr_MA[i]
        if ((operation_graph[i][0] <= tau_s <= operation_graph[i][1]) or (operation_graph[i][0] <= tau_s+tau_d <= operation_graph[i][1])) and operation_graph[i][2] == tau_m:
            operation_graph[i][0] = tau_s+tau_d
            operation_graph[i][1] = tau_s+tau_d+benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]]
            operation_graph[i][2] = opr_MA[i]
        machine_available_time[opr_MA[i]] = operation_graph[i][1]     #機械の操作可能時間の更新
        job_available_time[opr_OS[i]] = operation_graph[i][1]         #仕事の割当て可能時間の更新
        plt.barh(y=operation_graph[i][2],width=benchmark[opr_list[opr_OS[i]]+operaton_cnt[opr_OS[i]]][opr_MA[i]],left=operation_graph[i][0],color=cmap(opr_OS[i]),ec='k')    #プロット
        plt.text((operation_graph[i][0]+operation_graph[i][1])/2,operation_graph[i][2],opr_OS[i],ha="center")
        operaton_cnt[opr_OS[i]] += 1        #各仕事の完了した操作数を1増やす
    plt.barh(y=tau_m,width=tau_d,left=tau_s,color="yellow",ec="k")
    plt.xlabel("Makespan")
    plt.ylabel("machine")
    plt.show()          #図の提示
    return max(job_available_time)

#ベンチマーク問題を読み込み、工程表を作る
f = open('MK2.txt','r',encoding='utf-8')    #ファイル読み込み
# print(f)
mat = []                                    #ベンチマーク問題
with f as fin:                              #ベンチマーク問題の読み込み
    for line in fin.readlines():
        row = []
        toks = line.split(' ')
        # print(toks)
        for tok in toks:
            num = int(tok)
            row.append(num)
        mat.append(row)
N_job = mat[0][0]                   #仕事数
N_machine = mat[0][1]               #機械数
num_opr = []                        #各仕事の操作数
for i in range(1,N_job+1):          #各仕事の操作数を保存
    num_opr.append(mat[i][0])
sum_opr = 0                         #全操作数を保存する変数
for i in range(N_job):
    sum_opr += num_opr[i]
mat_ope = [[10**5]*N_machine for i in range(sum_opr)]     #各機械の操作時間表
sum_opr_i = 0                                         #表にまとめた操作数
opr_list = []                                         #工程表を仕事ごとに区切るためのリスト
for i in range(N_job):                                #工程表をまとめる
    f = 1                                             #使える機械数が書いてある列
    available_machine_num = mat[i+1][f]
    for j in range(num_opr[i]):
        for k in range(available_machine_num):      
            mat_ope[sum_opr_i+j][mat[i+1][f+2*k+1]-1] = mat[i+1][f+2*k+2]   #工程表に時間を書き込み 
        f = f+2*available_machine_num+1             #列の更新
        if j<num_opr[i]-1:                          #最終の操作でなければ次の操作で使える機械数を更新
            available_machine_num = mat[i+1][f]
    opr_list.append(sum_opr_i)
    sum_opr_i += num_opr[i]                         #工程表にまとめた機械数を更新
# pprint.pprint(mat_ope)

# for h in range(10):
#各データの空リスト
population_OS = []                      #操作順序セット
population_MA = []                      #割当て機械セット
population_size = 5*N_job*N_machine     #集団サイス
sol_MA = []
sol_OS = []
elite_ma = []
elite_os = []
list_cross = list(range(population_size))       #遺伝子のペア
min_list = []                                   #世代ごとのメイクスパンの最小値を格納するリスト
min_list_improvement = []                       #世代ごとのメイクスパンのずれの最小値を格納するリスト
t_list = []
#初期集団の生成
for i in range(population_size):
    sol_OS, sol_MA = initialization(num_opr,N_machine,opr_list,mat_ope,N_job)     #初期化
    population_OS.append(sol_OS)                                    #解の代入
    population_MA.append(sol_MA)                                    #解の代入
makespan_gene = [0 for i in range(population_size)]                 #メイクスパンを格納するリスト
evaluation_gene = [0 for i in range(population_size)]               #評価値を格納するリスト
makespan_gene_after_machine_break = [0 for i in range(population_size)]     #故障後のメイクスパンを格納するリスト
improvement_gene = [0 for i in range(population_size)]              #故障前と故障後のメイクスパンのずれを格納するリスト

#初期集団各個体のメイクスパンを計算
for j in range(population_size):                        
        makespan_gene[j] = makespan(population_MA[j],population_OS[j],N_machine,sum_opr,mat_ope,N_job,opr_list)
# print(min(makespan_gene))   #最小メイクスパンを出力

#パラメータ
p_c = 0.8                                       #交叉確率
p_m = 0.05                                      #突然変異確率
ganma = 0.6                                     #第二段階のメイクスパンとずれの割合
max_iteration = 10*N_job*N_machine              #最大反復回数
# max_iteration = 100
elite_os,elite_ma = elitizm(population_OS,population_MA,makespan_gene,population_size)                  #エリート個体の保存
#GA
#第一段階
for i in range(max_iteration):
    for j in range(population_size):            #評価値計算
        evaluation_gene[j] = 1-makespan_gene[j]/sum(makespan_gene)
    population_OS,population_MA = selection(population_OS,population_MA,evaluation_gene,population_size)    #選択
    r.shuffle(list_cross)                       #遺伝子のペアを作成
    for j in range(0,population_size,2):
        if r.random()<p_c:                      #交叉
            population_OS[list_cross[j]],population_OS[list_cross[j+1]],population_MA[list_cross[j]],population_MA[list_cross[j+1]] = pox(population_OS[list_cross[j]],population_OS[list_cross[j+1]],population_MA[list_cross[j]],population_MA[list_cross[j+1]],N_job,N_machine,sum_opr)
    for j in range(population_size):
        if r.random()<p_m:                      #突然変異(操作)
            population_OS[j],population_MA[j] = swap_mutation_os(population_OS[j],population_MA[j],len(population_OS[j]))
    # population_MA = mutation_ma(population_MA,population_OS,N_machine,num_opr,opr_list,p_m,mat_ope)
            population_MA[j] = mutation_ma(population_MA[j],population_OS[j],N_machine,num_opr,opr_list,mat_ope)     #突然変異(機械)
    for j in range(population_size):            #メイクスパンの評価
        makespan_gene[j] = makespan(population_MA[j],population_OS[j],N_machine,sum_opr,mat_ope,N_job,opr_list)
    for j in range(int(population_size/10)):    #エリート保存戦略
        s = makespan_gene.index(max(makespan_gene))     #メイクスパン最大の個体のインデックスを出力
        # if makespan_gene[s] > makespan(elite_ma[j],elite_os[j],N_machine,sum_opr,mat_ope,N_job,opr_list):
        population_OS[s] = elite_os[j]                  #エリート個体を代入
        population_MA[s] = elite_ma[j]
        makespan_gene[s] = makespan(population_MA[s],population_OS[s],N_machine,sum_opr,mat_ope,N_job,opr_list)
    # print(min(makespan_gene))
    min_list.append(min(makespan_gene))     #最良値を保存
    t_list.append(i)                        #世代数を保存
    elite_os,elite_ma = elitizm(population_OS,population_MA,makespan_gene,population_size)  #エリート主義
# print(min(makespan_gene))       #最小のメイクスパンを出力
# plt.plot(t_list,min_list)       #解の推移を描画
# plt.xlabel("generation")
# plt.ylabel("makespan")
# plt.show()                      #描画の出力
# g_d = makespan_gene.index(min(makespan_gene))     #最小メイクスパンの個体のインデックスを求める
# ganttchart(population_MA[g_d],population_OS[g_d],N_machine,sum_opr,mat_ope,N_job,opr_list)      #ガントチャートの作成
g = makespan_gene.index(min(makespan_gene))
a = r.uniform(0,0.5)
b = r.uniform(0.1,0.15)
DMS, DSTB = calculation_MS_and_STB(population_MA[g],population_OS[g],N_machine,sum_opr,mat_ope,N_job,opr_list,a,b)
print(DMS)
print(DSTB)

Z_gene = [0 for i in range(population_size)]
#初期集団の生成(第二段階)
for j in range(population_size):                     
    a = r.uniform(0,0.5)
    b = r.uniform(0.1,0.15)
    MS, STB = calculation_MS_and_STB(population_MA[j],population_OS[j],N_machine,sum_opr,mat_ope,N_job,opr_list,a,b)
    Z_gene[j] = ganma*MS+(1-ganma)*STB
elite_os,elite_ma = elitizm(population_OS,population_MA,Z_gene,population_size)                  #エリート個体の保存

# max_iteration = 100     #第二段階の世代数
max_iteration = 10*N_job*N_machine              #最大反復回数
#第二段階
for i in range(max_iteration):
    for j in range(population_size):            #評価値計算
        evaluation_gene[j] = 1-Z_gene[j]/sum(Z_gene)
    population_OS,population_MA = selection(population_OS,population_MA,evaluation_gene,population_size)    #選択
    r.shuffle(list_cross)                       #遺伝子のペアを作成
    for j in range(0,population_size,2):
        if r.random()<p_c:                      #交叉
            population_OS[list_cross[j]],population_OS[list_cross[j+1]],population_MA[list_cross[j]],population_MA[list_cross[j+1]] = pox(population_OS[list_cross[j]],population_OS[list_cross[j+1]],population_MA[list_cross[j]],population_MA[list_cross[j+1]],N_job,N_machine,sum_opr)
    for j in range(population_size):
        if r.random()<p_m:                      #突然変異(操作)
            population_OS[j],population_MA[j] = swap_mutation_os(population_OS[j],population_MA[j],len(population_OS[j]))
    # population_MA = mutation_ma(population_MA,population_OS,N_machine,num_opr,opr_list,p_m,mat_ope)
            population_MA[j] = mutation_ma(population_MA[j],population_OS[j],N_machine,num_opr,opr_list,mat_ope)     #突然変異(機械)
    for j in range(population_size):            #ずれの評価
        a = r.uniform(0,0.5)
        b = r.uniform(0.1,0.15)
        MS, STB = calculation_MS_and_STB(population_MA[j],population_OS[j],N_machine,sum_opr,mat_ope,N_job,opr_list,a,b)
        Z_gene[j] = ganma*MS+(1-ganma)*STB
    for j in range(int(population_size/10)):    #エリート保存戦略
        s = Z_gene.index(max(Z_gene))     #目的関数が最大の個体のインデックスを出力
        # if makespan_gene[s] > makespan(elite_ma[j],elite_os[j],N_machine,sum_opr,mat_ope,N_job,opr_list):
        population_OS[s] = elite_os[j]                  #エリート個体を代入
        population_MA[s] = elite_ma[j]
        # improvement_gene[s] = makespan(population_MA[s],population_OS[s],N_machine,sum_opr,mat_ope,N_job,opr_list)
        a = r.uniform(0,0.5)
        b = r.uniform(0.1,0.15)
        MS, STB = calculation_MS_and_STB(population_MA[s],population_OS[s],N_machine,sum_opr,mat_ope,N_job,opr_list,a,b)
        Z_gene[s] = ganma*MS+(1-ganma)*STB
    # print(min(makespan_gene))
    min_list.append(min(Z_gene))
    # min_list_improvement.append(min(improvement_gene))     #最良値を保存
    t_list.append(i)                        #世代数を保存
    elite_os,elite_ma = elitizm(population_OS,population_MA,Z_gene,population_size)  #エリート主義
# print(min(Z_gene))       #最小のメイクスパンを出力
# plt.plot(t_list,min_list)       #解の推移を描画
# plt.plot(t_list,min_list_improvement)
# plt.xlabel("generation")
# plt.ylabel("Z")
# plt.show()                      #描画の出力
g = Z_gene.index(min(Z_gene))     #最小のずれの個体のインデックスを求める
# print(makespan(population_MA[g],population_OS[g],N_machine,sum_opr,mat_ope,N_job,opr_list)) #最小のずれの個体のメイクスパンを求める
# ganttchart(population_MA[g],population_OS[g],N_machine,sum_opr,mat_ope,N_job,opr_list)      #ガントチャートの作成(機械の故障なし)
#ガントチャートの作成(機械の故障あり)
a = r.uniform(0,0.5)
b = r.uniform(0.1,0.15)
RMS,RSTB = calculation_MS_and_STB(population_MA[g],population_OS[g],N_machine,sum_opr,mat_ope,N_job,opr_list,a,b)
print(RMS)
print(RSTB)
# tau_m, MBT = select_break_machine(population_OS[g],population_MA[g],N_machine,sum_opr,N_job,mat_ope,opr_list)
# tau_s = a*MBT
# tau_d = b*MBT
# RMS = ganttchart_machine_break(population_MA[g],population_OS[g],N_machine,sum_opr,mat_ope,N_job,opr_list,tau_m,tau_s,tau_d)
# print(RMS)
# RMS, RSTB = calculation_MS_and_STB(population_MA[g],population_OS[g],N_machine,sum_opr,mat_ope,N_job,opr_list,a,b)
# print(RMS)

# Z_gene[s] = ganma*MS+(1-ganma)*STB
