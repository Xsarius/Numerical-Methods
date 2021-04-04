using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace GdGC
{
    public partial class Form1 : Form
    {

        private GD_Relation[] relation;

        static string dataFile = @".\Data.txt";

        public Form1()
        {
            InitializeComponent();
            LoadData(dataFile);

            comboBox1.SelectedIndex = 0;
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }

        private void textBox1_TextChanged(object sender, EventArgs e)
        {

        }

        private void textBox2_TextChanged(object sender, EventArgs e)
        {
                
        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
           // label1.Text = relation[comboBox1.SelectedIndex].yLabel;
            //label2.Text = relation[comboBox1.SelectedIndex].xLabel;
        }

        private void LoadData(string filePath)
        {
            
        }
    }
}
