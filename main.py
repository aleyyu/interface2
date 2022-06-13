import sys
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import QDir
from PyQt5.QtWidgets import *
from pathlib import Path
import pyqtgraph as pg
from pyqtgraph import PlotWidget

from interface import Ui_MainWindow
import fm_index_algorithm
import kmer_algorithm
import aho_corasick_algorithm

class MainWindow:
    def __init__(self):
        self.main_window = QMainWindow()    #QMainWindow instance
        self.ui = Ui_MainWindow()           #Arayüzü kullanabilmek için instance
        self.ui.setupUi(self.main_window)

        self.ui.stackedWidget.setCurrentWidget(self.ui.main_page)  #Stacked widget'ın current widget'ı mainwindow yapıldı.

        #butonların fonksiyonları atandı
        self.ui.btn_main.clicked.connect(lambda : self.goToPage(self.ui.main_page))
        self.ui.btn_bloom_aho.clicked.connect(lambda : self.goToPage(self.ui.bloom_aho_page))
        self.ui.btn_fm_kmer.clicked.connect(lambda : self.goToPage(self.ui.fm_page))
        self.ui.btn_bloom_aho.clicked.connect(lambda : self.goToPage(self.ui.bloom_aho_page))
        self.ui.btn_smith.clicked.connect(lambda : self.goToPage(self.ui.smith_page))

        #FM INDEX BUTTONS
        self.ui.fm_btn_upload_txt.clicked.connect(self.fm_load_txt)
        self.ui.fm_btn_upload_ptrn.clicked.connect(self.fm_load_ptrn)
        self.ui.fm_btn_search.clicked.connect(self.fm_search)
        self.ui.fm_btn_clear.clicked.connect(self.clear)
        self.ui.fm_btn_calc_kmer.clicked.connect(self.calc_kmer)
        self.ui.fm_btn_calc_kmer.clicked.connect(self.btn_kmer)
        self.ui.fm_btn_calc_kmer.clicked.connect(self.kmer_freq_top)
        #self.ui.fm_btn_calc_kmer.clicked.connect(lambda : self.kmer_draw_graph_1())
        #self.ui.fm_btn_calc_kmer.clicked.connect(lambda: self.kmer_draw_graph_2())
        #self.ui.fm_btn_calc_kmer.clicked.connect(lambda: self.kmer_draw_graph_3())


        #AHO CORASICK BUTTONS
        self.ui.aho_btn_upload_txt.clicked.connect(self.aho_load_txt)
        self.ui.aho_btn_upload_ptrn.clicked.connect(self.aho_load_ptrn)
        self.ui.bloom_btn_upload_txt.clicked.connect(self.bloom_load_txt)
        self.ui.bloom_btn_upload_ptrn.clicked.connect(self.bloom_load_ptrn)

        #BLOOM FILTER BUTTONS

        #SMITH WATERMAN BUTTONS



    def show(self):
        self.main_window.show()

    def goToPage(self, widget):
        self.ui.stackedWidget.setCurrentWidget(widget)

    #FM INDEX İÇİN GEREKLİ KODLAR BÖLÜMÜ
    txt = ""
    ptrn = ""
    sonuc = []
    k = 0
    f_name_txt = ""
    f_name_ptrn = ""

    def clear(self):
        self.ui.fm_output.clear()
        self.ui.fm_kmer_output.clear()
        self.ui.fm_kmer_output_top5.clear()
        self.ui.fm_graph1.clear()
        self.ui.fm_graph2.clear()
        self.ui.fm_graph3.clear()

    # ALIGNMENT
    def fm_search(self):
        t = self.txt
        p = self.ptrn
        count = 0
        output = []
        run_time = 0.0

        if self.f_name_txt == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for text sequence selected. Please select one.")
            QMessageBox.show()

        elif self.f_name_ptrn == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for pattern sequence selected. Please select one.")
            QMessageBox.show()
        else:
            result = fm_index_algorithm.align(t, p)
            output = result[0]
            run_time = result[1]
            #print(result)

            self.ui.fm_output.setText("The start position of pattern sequence is: ")

            for i in output:
                # self.textBrowser_3.append("i: {}".format(i))
                self.ui.fm_output.append("- {}".format(i))
                count += 1

            self.ui.fm_output.append("\n-----------------------------------------------------------------------------------------------")
            self.ui.fm_output.append("\n {} start position found".format(count))
            self.ui.fm_output.append("\n-----------------------------------------------------------------------------------------------")
            self.ui.fm_output.append("\n Run time is: {}".format(run_time))


    def fm_load_txt(self):
        dialog = QFileDialog(directory = "D:\Masaüstü\FASTA files")
        dialog.setWindowTitle("Choose FASTA file for text sequence.")
        dialog.setNameFilter("FNA files (*.fna)")
        dialog.setFilter(QDir.Files)

        str_fname = ""

        if dialog.exec_():
            file_name = dialog.selectedFiles()
            self.f_name_txt = file_name

            if file_name[0].endswith('.fna'):
                with open(file_name[0], 'r') as f:
                    text = f.read()
                    self.txt = text
                    f.close()
            else:
                pass

        str_fname = str_fname.join(file_name)


        if file_name:
            fname = Path(str_fname)
            self.ui.fm_lbl_txt.setText("Selected file: {}".format(fname.name))
        else:
            pass
            self.ui.fm_lbl_txt.setText("No file is chosen. Please choose a file.")

            return text

    def fm_load_ptrn(self):
        dialog = QFileDialog(directory = "D:\Masaüstü\FASTA files")
        dialog.setWindowTitle("Choose FASTA file for pattern sequence.")
        dialog.setNameFilter("FNA files (*.fna)")
        dialog.setFilter(QDir.Files)

        str_fname = ""

        if dialog.exec_():
            file_name = dialog.selectedFiles()
            self.f_name_ptrn = file_name

            if file_name[0].endswith('.fna'):
                with open(file_name[0], 'r') as f:
                    pattern = f.read()
                    self.ptrn = pattern
                    f.close()
            else:
                pass

        str_fname = str_fname.join(file_name)

        if file_name:
            fname = Path(str_fname)
            self.ui.fm_lbl_ptrn.setText("Selected file: {}".format(fname.name))
        else:
            self.ui.fm_lbl_ptrn.setText("No file was choosed. Please choose a file.")

            return pattern

    #KMER
    def btn_kmer(self):
        sequence = self.ptrn
        get_k = self.ui.fm_input_kmer_value.text()

        if self.f_name_ptrn == "":
            QMessageBox.about(main_window, "ERROR", "No FASTA file for pattern sequence selected. Please select one.")
            QMessageBox.show()

        elif get_k == "":
            QMessageBox.about(main_window, "ERROR", "K value for K-Mer algorithm can't be 0. Please enter a valid value.")
            QMessageBox.show()
        else:
            self.calc_kmer()
            self.kmer_freq_top()
            self.draw_graph_1()
            self.draw_graph_2()
            self.draw_graph_3()

    def calc_kmer(self):
        sequence = self.ptrn
        get_k = self.ui.fm_input_kmer_value.text()

        ksize = int(get_k)
        self.k = ksize
        kmers = kmer_algorithm.build_kmers(sequence, ksize)

        for kmer, freq in kmers.items():
            txt = "kmer: {} - freq: {}".format(kmer, freq)
            self.ui.fm_kmer_output.append(txt)
        #verScrollBar = self.txt_kmer_freq.setVerticalScrollBar()
        #verScrollBar.setValue(verScrollBar.minimum())

    def kmer_freq_top(self):

        sequence = self.ptrn
        get_k = self.ui.fm_input_kmer_value.text()
        ksize = int(get_k)
        self.k = ksize

        mc = kmer_algorithm.kmer_freq_top(sequence, ksize)

        for kmer, freq in mc:
            txt = "kmer: {} - freq: {}".format(kmer, freq)
            self.ui.fm_kmer_output_top5.append("{}".format(txt))
            print("{}".format(txt))

    def kmer_draw_graph_1(self):
        sequence = self.ptrn
        get_k = self.ui.fm_input_kmer_value.text()
        ksize = int(get_k)
        self.k = ksize

        x = []
        y = []

        sonuc = kmer_algorithm.kmer_graph_1(sequence, ksize)
        x = sonuc[0]
        y = sonuc[1]

        pen = pg.mkPen(color=(255, 0, 0), width=1, style=QtCore.Qt.DotLine)
        styles = {'color': 'r', 'font-size': '15px'}
        self.ui.fm_graph1.setLabel('left', '1/freq', **styles)
        self.ui.fm_graph1.setLabel('bottom', 'freqs', **styles)
        self.ui.fm_graph1.setTitle("K-Mer Period Graph",color="red", size="10pt")
        self.ui.fm_graph1.plot(x, y, pen=pen)
        self.ui.fm_graph1.updateMatrix()

    def kmer_draw_graph_2(self):
        sequence = self.ptrn
        get_k = self.ui.fm_input_kmer_value.text()
        ksize = int(get_k)
        self.k = ksize

        x = []
        y = []

        sonuc2 = kmer_algorithm.kmer_graph_2(sequence, ksize)
        x = sonuc2[0]
        y = sonuc2[1]

        #pen = pg.mkPen(color="blue", width=1, style=QtCore.Qt.SolidLine)
        styles = {'color': 'r', 'font-size': '15px'}
        self.ui.fm_graph2.setLabel('left', 'freq', **styles)
        self.ui.fm_graph2.setLabel('bottom', 'k-mer no', **styles)
        self.ui.fm_graph2.setTitle("K-Mer Frequence Graph", color="b", size="10pt")
        #self.graph_2.plot(x, y,pen=pen, symbol='+', symbolSize=5, symbolBrush=('b'))
        self.ui.fm_graph2.plot(x, y, pen=None, symbol='o')
        self.ui.fm_graph2.updateMatrix()


    def kmer_draw_graph_3(self):
        sequence = self.ptrn
        get_k = self.ui.fm_input_kmer_value.text()
        ksize = int(get_k)
        self.k = ksize

        x = []
        y = []

        sonuc3 = kmer_algorithm.kmer_graph_3(sequence, ksize)
        x = sonuc3[0]
        y = sonuc3[1]

        #pen = pg.mkPen(color="white", width=1, style=QtCore.Qt.SolidLine)
        bargraph = pg.BarGraphItem(x=x, height=y, width=0.6, brush='g')
        styles = {'color': 'r', 'font-size': '15px'}
        self.ui.fm_graph3.setLabel('left', 'freq', **styles)
        self.ui.fm_graph3.setLabel('bottom', 'kmer no', **styles)
        self.ui.fm_graph3.setTitle("Top 5 Frequence Graph", color="g", size="10pt")
        self.ui.fm_graph3.plot(x,y)
        self.ui.fm_graph3.addItem(bargraph)
        self.ui.fm_graph3.updateMatrix()

    a_ptrn_fnames = []

    #AHO CORASICK KODLARI
    def aho_load_txt(self):

        dialog = QFileDialog(directory = "D:\Masaüstü\FASTA files")
        dialog.setWindowTitle("Choose FASTA file for text sequence.")
        dialog.setNameFilter("FNA files (*.fna)")
        dialog.setFilter(QDir.Files)
        dialog.setFileMode(QFileDialog.ExistingFiles)

        str_fname = ""

        if dialog.exec_():
            file_name = dialog.selectedFiles()
            self.f_name_txt = file_name

            if file_name[0].endswith('.fna'):
                with open(file_name[0], 'r') as f:
                    text = f.read()
                    self.txt = text
                    f.close()
            else:
                pass

        str_fname = str_fname.join(file_name)


        if file_name:
            fname = Path(str_fname)
            self.ui.fm_lbl_txt.setText("Selected file: {}".format(fname.name))
        else:
            pass
            self.ui.fm_lbl_txt.setText("No file is chosen. Please choose a file.")

            return text

    def aho_load_ptrn(self):

        caption = "Choose FASTA file for pattern sequence."
        dir = 'D:\Masaüstü\FASTA files'
        filter = "FNA files (*.fna)"
        str_fname = ""

        file_names = QFileDialog.getOpenFileNames(None, caption, dir, filter)[0]
        for file in file_names:
            self.a_ptrn_fnames = file
            print(file)
            if file.endswith(".fna"):
                with open(file, 'r') as f:
                    pattern = f.read()
                    self.ptrn = pattern
                    f.close()
            else:
                pass

            """
            str_fname = str_fname.join(file)

            if file:
                fname = Path(str_fname)
                self.ui.aho_lbl_ptrn.setText("Selected file: {}".format(fname.name))
            else:
                self.ui.aho_lbl_ptrn.setText("No file was choosed. Please choose a file.")
                return pattern
            """


    def bloom_load_txt(self):
        pass

    def bloom_load_ptrn(self):
        pass


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())