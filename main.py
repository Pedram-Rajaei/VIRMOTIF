import os
from tkinter import *
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import axes3d
import numpy as np
from matplotlib import cm
import pandas
from sklearn.cluster import KMeans
# from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
from tkinter import filedialog
from sklearn.preprocessing import StandardScaler
from functools import cmp_to_key
from sklearn.cluster import MeanShift, estimate_bandwidth
from itertools import cycle
import seaborn as sns

LARGE_FONT= ("Verdana", 12)

class Data:
    _instance = None
    @staticmethod
    def getInstance():
        if Data._instance == None:
            Data()
        else:
            return Data._instance
    def __init__(self):
        if Data._instance != None:
            raise Exception("multiple Contstruction!")
        else:
            self.file = ''
            self.values = [{}]
            self.probs = [{}]
            self.allCode = []
            self.mode = None
            Data._instance = self
    def setFile(self, file):
        self.file = file
    def setV(self, v):
        self.v = v
    def makeStr(self):
        self.allCode = ['A', 'C', 'G', 'T']
        i = 0
        while True:
            if len(self.allCode[i]) == 5:
                break
            for c in ['A', 'C', 'G', 'T']:
                self.allCode.append(self.allCode[i] + c)
            i = i+1
    def makeFile(self):
        self.makeStr()
        with open('./result.csv', 'w') as f:
            f.write('Sequence Number')
            for code in self.allCode:
                f.write(',' + code)
            f.write('\n')
        with open('./results_normalize.csv', 'w') as f:
            f.write('Sequence Number')
            for code in self.allCode:
                f.write(',' + code)
            f.write('\n')
    def addSequence(self, s, index):
        self.values.append({})
        self.probs.append({})
        cnt = [{}]
        sums = [0 for i in range(6)]
        for j in range(1, 6):
            cnt.append({})
            for i in range(len(s)-j):
                ss = s[i:i+j]
                if len(ss) != j:
                    print('bad index: ', j)
                if cnt[j].get(ss) == None:
                    cnt[j][ss] = 0
                cnt[j][ss] += 1
                sums[j] += 1
        for code in self.allCode:
            if cnt[len(code)].get(code) == None:
                self.values[index][code] = 0
                self.probs[index][code] = 0
            else:
                self.values[index][code] = cnt[len(code)][code]
                self.probs[index][code] = cnt[len(code)][code] / sums[len(code)]
        with open('result.csv', 'a') as f:
            with open('results_normalize.csv', 'a') as g:
                f.write('Seq ' + str(index))
                g.write('Seq ' + str(index))
                for code in self.allCode:
                    f.write(',' + str(self.values[index][code]))
                    g.write(',' + str(self.probs[index][code]))
                f.write('\n')
                g.write('\n')
    def calculate(self):
        index = 0
        with open(self.file, 'r') as f:
            lines = f.readlines()
            s = ''
            for line in lines:
                line = line.strip()
                if(len(line)>0 and line[0] == '>'):
                    if len(s) > 0:
                        self.addSequence(s, index)
                    s = ''
                    index = index + 1
                    # if index > 5:
                    #     break
                else:
                    for l in line:
                        if l in ['A', 'C', 'G', 'T']:
                            s = s + l
    def getSubstrings(self, sub_len, subs):
        results = []
        for index, code in enumerate(self.allCode):
            if len(code) == sub_len or code in subs:
                results.append([index, code])
        return results
    def getFrequencies(self, code):
        results = []
        for i in range(1, len(self.values)):
            if self.values[i][code] != None:
                results.append([i, self.values[i][code]])
        return results
    def getValues(self, code):
        if self.v == '1':
            return self.getFrequencies(code)
        else:
            return self.getDeRatio(code)
    def getProb(self, index, code):
        if len(code) > 0:
            return self.probs[index][code]
        else:
            return 1
    def getDeRatio(self, code):
        results = []
        for i in range(1, len(self.values)):
            if len(code) == 1:
                results.append([i, self.getProb(i, code)])
            else:
                results.append([i, self.getProb(i, code[:-1]) * self.getProb(i, code[1:]) / self.getProb(i, code[1:-1])])
        return results
data = Data()

class Three:
    def __init__(self, s, x, y):
        self.s = s
        self.x = x
        self.y = y
    def __lt__(self, o):
        if self.x == o.x:
            return self.y < o.y
        return self.x < o.x

class App(Tk):
    def __init__(self, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)
        container = Frame(self)
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        self.sub_len = 0
        self.subs = []
        self.frames = {}
        self.mode = None
        for F in (StartPage, PagePlot, PageClustering, SelectLengthPage, SelectSubstringPage):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")
        self.show_frame(StartPage)
    def show_frame(self, cont):
        self.frames[cont].update()
        frame = self.frames[cont]
        frame.tkraise()
    def show(self, seqs, subs):
        if plt != None:
            plt.close()
        names = []
        for s in seqs:
            s = 'Seq ' + s
            names.append(s)

        result = []
        for seq in seqs:
            a = []
            for sub in subs:
                if data.v == '1':
                    a.append(str(data.values[int(seq)][sub]))
                else:
                    a.append(str(data.probs[int(seq)][sub]))
            result.append(a)

        result = np.array(result, dtype=np.float)

        fig=plt.figure(figsize=(5, 5), dpi=150)
        ax1=fig.add_subplot(111, projection='3d')

        xlabels = np.array(subs)
        ylabels = np.array(names)
        xpos = np.arange(xlabels.shape[0])
        ypos = np.arange(ylabels.shape[0])

        xposM, yposM = np.meshgrid(xpos, ypos, copy=False)

        zpos=result
        zpos = zpos.ravel()

        dx=0.5
        dy=0.5
        dz=zpos

        ax1.w_xaxis.set_ticks(xpos + dx/2.)
        ax1.w_xaxis.set_ticklabels(xlabels)

        ax1.w_yaxis.set_ticks(ypos + dy/2.)
        ax1.w_yaxis.set_ticklabels(ylabels)

        values = np.linspace(0.2, 1., xposM.ravel().shape[0])
        colors = cm.rainbow(values)
        ax1.bar3d(xposM.ravel(), yposM.ravel(), dz*0, dx, dy, dz, color=colors)
        plt.show()
    def meanshift(self, subs, sub_len):
        a = []
        y = []
        for index, code in data.getSubstrings(sub_len, subs):
            for i, value in data.getValues(code):
                a.append([index, value])
                y.append(code)
        X = np.array(a)
        bandwidth = estimate_bandwidth(X, quantile=0.2, n_samples=500)
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(X)
        labels = ms.labels_
        cluster_centers = ms.cluster_centers_
        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)
        plt.figure(1)
        plt.clf()
        colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
        for k, col in zip(range(n_clusters_), colors):
            my_members = labels == k
            cluster_center = cluster_centers[k]
            plt.plot(X[my_members, 0], X[my_members, 1], col + '.')
            # plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
            #         markeredgecolor='k', markersize=14)
        plt.title('Estimated number of clusters: %d' % n_clusters_)
        plotted = []
        for i in range(X.shape[0]):
            if y[i] not in plotted:
                plt.annotate(y[i], (float(X[i, 0]), float(X[i, 1])))
                plotted.append(y[i])
        plt.show()
    def kmeans(self, subs, sub_len, k):
        a = []
        y = []
        for index, code in data.getSubstrings(sub_len, subs):
            for i, value in data.getValues(code):
                a.append([index, value])
                y.append(code)
        X = np.array(a)
        kmeans = KMeans(n_clusters=k, random_state=0).fit(X)
        labels = kmeans.labels_
        if plt != None:
            plt.close()
        plt.scatter(X[:, 0], X[:, 1], c = labels)
        plotted = []
        for i in range(X.shape[0]):
            if y[i] not in plotted:
                plt.annotate(y[i], (float(X[i, 0]), float(X[i, 1])))
                plotted.append(y[i])
        plt.show()
    def heatmap(self, subs, sub_len):
        a = []
        y = []
        for index, code in data.getSubstrings(sub_len, subs):
            for _, value in data.getValues(code):
                a.append([index, value])
                y.append(code)
        if plt != None:
            plt.close()
        ax1 = plt.axes()
        X = np.array(a)
        sns.set(color_codes=True)
        g = sns.heatmap(X, ax=ax1)
        plt.show()
    def clustermap(self, subs, sub_len):
        a = []
        y = []
        for index, code in data.getSubstrings(sub_len, subs):
            for _, value in data.getValues(code):
                a.append([index, value])
                y.append(code)
        if plt != None:
            plt.close()
        X = np.array(a)
        sns.set(color_codes=True)
        g = sns.clustermap(X, figsize=(5, 5))
        plt.show()
    def cmpByX(self, x, y):
        if x[0] != y[0]:
            return x[0] - y[0]
        return x[1] - y[1]
    def findIndex(self, sortedX, x):
        low, high = 0, len(sortedX)
        while (low+1) < high:
            mid = (low + high) // 2
            if sortedX[mid].x > x:
                high = mid
            else:
                low = mid
        return low
    def findNearest(self, sortedX, sortedY, x, y, T, debug=False):
        low = self.findIndex(sortedX, x)
        st = self.findIndex(sortedX, x-T)
        en = self.findIndex(sortedX, x+T)
        ans = 0
        for i in range(st, en):
            if (abs(x - sortedX[i].x) + abs(y - sortedX[i].y)) < T:
                ans = ans + 1
        return ans
    def pca(self, subs, sub_len, k):
        a = []
        L = []
        for index, code in data.getSubstrings(sub_len, subs):
            for i, value in data.getValues(code):
                a.append([index, i, value])
                L.append(code)
        X = np.array(a)
        pca = PCA(n_components=2)
        X = pca.fit_transform(X)
        scaler = StandardScaler()
        scaler.fit(X)
        X=scaler.transform(X)    
        pca = PCA()
        x_new = pca.fit_transform(X)

        def myplot(score,coeff,labels=None):
            xs = score[:,0]
            ys = score[:,1]
            n = coeff.shape[0]
            scalex = 1.0/(xs.max() - xs.min())
            scaley = 1.0/(ys.max() - ys.min())
            kmeans = KMeans(n_clusters=k).fit(score)
            plt.scatter(xs * scalex, ys * scaley, c=kmeans.labels_)
            X = np.array([xs * scalex, ys * scaley])
            xs = xs * scalex
            ys = ys * scaley
            oX = []
            for x, y in zip(xs, ys):
                oX.append([x, y])
            X = np.around(X, decimals=3)
            oX = []
            for i in range(0, len(X[0, :])):
                oX.append(Three(L[i], float(X[0, i]), float(X[1, i])))
            X = np.array([L, X[0, :], X[1, :]]).T
            mn = 10**10
            mx = -10**10
            for x in X:
                mx = max(mx, float(x[1]))
                mn = min(mn, float(x[1]))
            T = .1 * (mx - mn)
            oX.sort()
            sortedX = oX
            for x in sortedX:
                number = self.findNearest(sortedX, sortedX, float(x.x), float(x.y), T, False)
                if number < 10:
                    plt.annotate(x.s, (float(x.x), float(x.y)))
            plt.xlim(-1,1)
            plt.ylim(-1,1)
            plt.xlabel("PC{}".format(1))
            plt.ylabel("PC{}".format(2))
            plt.grid()

        #Call the function. Use only the 2 PCs.
        myplot(x_new[:,0:2],np.transpose(pca.components_[0:2, :]))
        plt.show()
    def load(self):
        file = filedialog.askopenfilename(initialdir = "./", title = "Select file", filetypes = [("all files", "*.*")])
        data.setFile(file)
        data.makeFile()
        data.calculate()

class StartPage(Frame):
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        self.parent = parent
        self.loaded = False
        self.controller = controller
        self.label1 = Label(self, text='', width=30)
        self.label1.grid(row=0, column=0)
        Label(self, text='FileName: ', width=15).grid(row=1, column=0)
        Label(self, text='Actions:', width=15).grid(row=2, column=0)
        self.filelabel = Label(self, text='', width=40)
        self.filelabel.grid(row=1, column=1)
        b_load = Button(self, width=10)
        b_load['text'] = 'Load'
        b_load['command'] = lambda: self.load()
        b_load.grid(row=1, column=2)

        self.v = StringVar(self, "1")
        values = {
            "Frequency Based": "1",
            "DeRatio Based": "2"
        }
        i = 0
        for (text, value) in values.items():
            r = Radiobutton(self, text=text, variable=self.v, value=value)
            r.grid(row=2, column=i)
            i = i+1
        b_plot = Button(self, width=10)
        b_plot['text'] = 'Plot'
        b_plot['command'] = self.plot
        b_plot.grid(row=2, column=2)

        b_clustering = Button(self, width=10)
        b_clustering['text'] = 'Clustering'
        b_clustering['command'] = self.clustering
        b_clustering.grid(row=3, column=2)
    def plot(self):
        if self.loaded == False:
            self.label1['text'] = 'file must be loaded before action'
            return
        data.setV(self.v.get())
        self.controller.show_frame(PagePlot)
    def clustering(self):
        if self.loaded == False:
            self.label1['text'] = 'file must be loaded before action'
            return
        data.setV(self.v.get())
        self.controller.show_frame(PageClustering)
    def load(self):
        self.label1['text'] = ''
        self.filelabel['text'] = 'Loading...'
        self.controller.load()
        self.filelabel['text'] = data.file
        if len(self.filelabel['text']) > 25:
            self.filelabel['text'] = '...' + self.filelabel['text'][-20:]
        self.loaded = True
    def update(self):
        pass

class PagePlot(Frame):
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        self.parent = parent
        self.controller = controller
        self.label1 = Label(self, text='')
        self.label1.grid(row=1, column=1)

        Label(self, text='Sequence Number', width=20).grid(row=2, column=0)
        self.e1 = Entry(self)
        self.e1.grid(row=2, column=1)

        Label(self, text='SubStrings', width=20).grid(row=3, column=0)
        self.e2 = Entry(self)
        self.e2.grid(row=3, column=1)

        Label(self, text='SubStrings Length', width=20).grid(row=4, column=0)
        self.e3 = Entry(self)
        self.e3.grid(row=4, column=1)
        b4 = Button(self, width=2)
        b4['text'] = '+'
        b4['command'] = self.selectLength
        b4.grid(row=4, column=2)


        Label(self, text='').grid(row=5, column=0)

        b1 = Button(self)
        b1['text'] = 'Bar chart'
        b1.grid(row=6, column=2)
        b1['command'] = lambda: self.plot()
        b2 = Button(self, width=10)
        b2['text'] = 'HeatMap'
        b2.grid(row=6, column=3)
        b2['command'] = lambda: self.heatmap()

        b3 = Button(self)
        b3['text'] = 'Back'
        b3.grid(row=6, column=0)
        b3['command'] = lambda: controller.show_frame(StartPage)
    def plot(self):
        if len(self.e1.get()) > 0:
            seqs = self.e1.get().split(',')
            seqs = [s.strip() for s in seqs]
        else:
            self.label1['text'] = 'no sequence number'
            return
        if len(self.e2.get()) > 0:
            subs = self.e2.get().split(',')
            subs = [s.strip() for s in subs]
        else:
            subs = []
        if len(self.e3.get()) > 0:
            subsE = [int(s.strip()) for s in self.e3.get().split(',')]
            for x in subsE:
                if x < 1 or x > 5:
                    self.label1['text'] = 'too short or too long substring'
                    return
            for sub in data.allCode:
                if len(sub) in subsE:
                    subs.append(sub)
        for seq in seqs:
            if len(data.values) <= int(seq):
                self.label1['text'] = 'invalid sequence number'
                return
        if len(subs) == 0:
            self.label1['text'] = 'empty substring list'
            return
        for sub in subs:
            for x in sub:
                if x not in ['A', 'C', 'G', 'T']:
                    self.label1['text'] = 'invalud character in substring'
                    return
        self.controller.show(seqs, subs)
    def heatmap(self):
        subs = []
        sub_len = 0
        if len(self.e2.get()) > 0:
            subs = self.e2.get().split(',')
            subs = [s.strip() for s in subs]
        if len(self.e3.get()) > 0:
            sub_len = int(self.e3.get())
        if len(self.e1.get()) > 0:
            self.label1['text'] = 'bad input: sequence number'
            return
        if sub_len == 0 and len(subs) == 0:
            self.label1['text'] = 'bad input: no substring'
            return
        self.controller.heatmap(subs, sub_len)
    def selectLength(self):
        self.controller.show_frame(SelectLengthPage)
    def update(self):
        self.e3.delete(0, END)
        self.e3.insert(END, str(self.controller.sub_len))
        self.controller.mode = 'Plot'

class PageClustering(Frame):
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        self.parent = parent
        self.controller = controller

        self.label1 = Label(self, text='')
        self.label1.grid(row=1, column=1)
        
        Label(self, text='').grid(row=2, column=1)

        Label(self, text='SubStrings', width=25).grid(row=3, column=0)
        self.e1 = Entry(self, width=25)
        self.e1.grid(row=3, column=1)

        Label(self, text='SubStrings Length', width=25).grid(row=4, column=0)
        self.e2 = Entry(self, width=25)
        self.e2.grid(row=4, column=1)
        b2 = Button(self, width=2)
        b2['text'] = '+'
        b2['command'] = self.selectLength
        b2.grid(row=4, column=2)

        Label(self, text='n-cluster/n-component', width=25).grid(row=5, column=0)
        self.e3 = Entry(self, width=25)
        self.e3.grid(row=5, column=1)

        Label(self, text='').grid(row=6, column=0)

        b3 = Button(self, width=10)
        b3['text'] = 'Back'
        b3['command'] = lambda: self.controller.show_frame(StartPage)
        b3.grid(row=7, column=0)
        b4 = Button(self, width=10)
        b4['text'] = 'K-means'
        b4.grid(row=7, column=2)
        b4['command'] = lambda: self.kmeans()
        b5 = Button(self, width=10)
        b5['text'] = 'PCA'
        b5.grid(row=7, column=3)
        b5['command'] = lambda: self.pca()
        b6 = Button(self, width=10)
        b6['text'] = 'MeanShift'
        b6.grid(row=8, column=2)
        b6['command'] = lambda: self.meanshift()
        b7 = Button(self, width=10)
        b7['text'] = 'ClusterMap'
        b7.grid(row=8, column=3)
        b7['command'] = lambda: self.clustermap()
    def selectSubstring(self):
        self.controller.show_frame(SelectSubstringPage)
    def selectLength(self):
        self.controller.show_frame(SelectLengthPage)
    def kmeans(self):
        subs = []
        sub_len = 0
        if len(self.e1.get()) > 0:
            subs = self.e1.get().split(',')
            subs = [s.strip() for s in subs]
        if len(self.e2.get()) > 0:
            sub_len = int(self.e2.get())
        if sub_len == 0 and len(subs) == 0:
            self.label1['text'] = 'bad input: no substring'
            return
        if len(self.e3.get()) > 0:
            k = int(self.e3.get())
        else:
            self.label1['text'] = 'n-cluster is null'
            return
        self.controller.kmeans(subs, sub_len, k)
    def clustermap(self):
        subs = []
        sub_len = 0
        if len(self.e1.get()) > 0:
            subs = self.e1.get().split(',')
            subs = [s.strip() for s in subs]
        if len(self.e2.get()) > 0:
            sub_len = int(self.e2.get())
        if sub_len == 0 and len(subs) == 0:
            self.label1['text'] = 'bad input: no substring'
            return
        self.controller.clustermap(subs, sub_len)
    def pca(self):
        subs = []
        sub_len = 0
        if len(self.e1.get()) > 0:
            subs = self.e1.get().split(',')
            subs = [s.strip() for s in subs]
        if len(self.e2.get()) > 0:
            sub_len = int(self.e2.get())
        if sub_len == 0 and len(subs) == 0:
            self.label1['text'] = 'bad input: no substring'
            return
        if len(self.e3.get()) > 0:
            k = int(self.e3.get())
        else:
            self.label1['text'] = 'n-cluster is null'
            return
        self.controller.pca(subs, sub_len, k)
    def meanshift(self):
        subs = []
        sub_len = 0
        if len(self.e1.get()) > 0:
            subs = self.e1.get().split(',')
            subs = [s.strip() for s in subs]
        if len(self.e2.get()) > 0:
            sub_len = int(self.e2.get())
        if sub_len == 0 and len(subs) == 0:
            self.label1['text'] = 'bad input: no substring'
            return
        self.controller.meanshift(subs, sub_len)
    def update(self):
        self.e2.delete(0, END)
        self.e2.insert(END, str(self.controller.sub_len))
        self.e1.delete(0, END)
        self.e1.insert(END, str(", ".join(self.controller.subs)))
        self.controller.mode = 'Clustering'

class SelectLengthPage(Frame):
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        self.parent = parent
        self.controller = controller

        self.label1 = Label(self, text='')
        self.label1.grid(row=0, column=0)

        self.vars = []
        for i in range(5):
            var = IntVar()
            Checkbutton(self, text=str(i+1), variable=var, width=10).grid(row=1, column=i)
            self.vars.append(var)

        Label(self, text='').grid(row=2, column=0)

        b1 = Button(self, width=10)
        b1['text'] = 'Add'
        b1['command'] = lambda: self.submit()
        b1.grid(row=3, column=2)
    def submit(self):
        cnt = 0
        idx = 0
        for i, v in enumerate(self.vars):
            cnt = cnt + int(v.get())
            if int(v.get()) == 1:
                idx = i
        if cnt != 1:
            self.label1['text'] = 'too many choices'
        else:
            self.controller.sub_len = idx+1
            if self.controller.mode == 'Plot':
                self.controller.show_frame(PagePlot)
            else:
                self.controller.show_frame(PageClustering)
    def update(self):
        pass

class SelectSubstringPage(Frame):
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        self.parent = parent
        self.controller = controller

        self.label1 = Label(self, text='')
        self.label1.grid(row=0, column=0)

        self.vars = []
        self.c = []
        for i, code in enumerate(data.allCode):
            var = IntVar()
            self.c.append(Checkbutton(self, text=code, variable=var, width=10).grid(row=i//10, column=i%10))
            self.vars.append(var)

        Label(self, text='').grid(row=2, column=0)

        b1 = Button(self, width=10)
        b1['text'] = 'Add'
        b1['command'] = lambda: self.submit()
        b1.grid(row=3, column=2)
    def submit(self):
        for i, v in enumerate(self.vars):
            if int(v.get()) == 1:
                self.controller.subs.append(self.c[i]['text'])
                self.controller.show_frame(PageClustering)
    def update(self):
        self.__init__(self.parent, self.controller)

app = App()
app.mainloop()
