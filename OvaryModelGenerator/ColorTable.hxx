#include <iostream>
#include <string>
#include <fstream>
#include <exception>

// The Color list from SNAP
// https://github.com/pyushkevich/itksnap/blob/312323ea9477ae8ee2c47dd20a93b621527d71d8/Logic/Common/ColorLabelTable.cxx#L39

const size_t ColorListSize = 130;
const char *ColorList[ColorListSize] = {
"#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#00FFFF", "#FF00FF",
"#FFEFD5", "#0000CD", "#CD853F", "#D2B48C", "#66CDAA", "#000080",
"#008B8B", "#2E8B57", "#FFE4E1", "#6A5ACD", "#DDA0DD", "#E9967A",
"#A52A2A", "#FFFAFA", "#9370DB", "#DA70D6", "#4B0082", "#FFB6C1",
"#3CB371", "#FFEBCD", "#FFE4C4", "#DAA520", "#008080", "#BC8F8F",
"#FF69B4", "#FFDAB9", "#DEB887", "#7FFF00", "#8B4513", "#7CFC00",
"#FFFFE0", "#4682B4", "#006400", "#EE82EE", "#EEE8AA", "#F0FFF0",
"#F5DEB3", "#B8860B", "#20B2AA", "#FF1493", "#191970", "#708090",
"#228B22", "#F8F8FF", "#F5FFFA", "#FFA07A", "#90EE90", "#ADFF2F",
"#4169E1", "#FF6347", "#FAF0E6", "#800000", "#32CD32", "#F4A460",
"#FFFFF0", "#7B68EE", "#FFA500", "#ADD8E6", "#FFC0CB", "#7FFFD4",
"#FF8C00", "#8FBC8F", "#DC143C", "#FDF5E6", "#FFFAF0", "#00CED1",
"#00FF7F", "#800080", "#FFFACD", "#FA8072", "#9400D3", "#B22222",
"#FF7F50", "#87CEEB", "#6495ED", "#F0E68C", "#FAEBD7", "#FFF5EE",
"#6B8E23", "#87CEFA", "#00008B", "#8B008B", "#F5F5DC", "#BA55D3",
"#FFE4B5", "#FFDEAD", "#00BFFF", "#D2691E", "#FFF8DC", "#2F4F4F",
"#483D8B", "#AFEEEE", "#808000", "#B0E0E6", "#FFF0F5", "#8B0000",
"#F0FFFF", "#FFD700", "#D8BFD8", "#778899", "#DB7093", "#48D1CC",
"#FF00FF", "#C71585", "#9ACD32", "#BDB76B", "#F0F8FF", "#E6E6FA",
"#00FA9A", "#556B2F", "#40E0D0", "#9932CC", "#CD5C5C", "#FAFAD2",
"#5F9EA0", "#008000", "#FF4500", "#E0FFFF", "#B0C4DE", "#8A2BE2",
"#1E90FF", "#F08080", "#98FB98", "#A0522D"};

struct LabelColor
{
  LabelColor() {};

  LabelColor(uint8_t _r, uint8_t _g, uint8_t _b, float _a)
    : r(_r), g(_g), b(_b), a(_a) {}

  LabelColor(const LabelColor& other)
  {
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
    this->a = other.a;
  }

  LabelColor& operator=(const LabelColor& other)
  {
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
    this->a = other.a;
    return *this;
  }

  uint8_t r, g, b;
  float a;
};


class ColorTable
{
public:
  typedef std::map<uint16_t, LabelColor> LabelColorMap;

  ColorTable() {} // default constructor using default colors
  ColorTable(const char* fnLabelDesc) // build color table using a label description file
  {
    // Use default color if file not provided
    if (strlen(fnLabelDesc) == 0)
      return; 

    // Create a stream for reading the file
    std::ifstream fin(fnLabelDesc);
    std::string line;

    // Create a temporary map of labels (to discard in case there is a problem
    // reading the file later)
    LabelColorMap inputMap;

    // Set the clear label in the input map
    //inputMap[0] = this->GetDefaultLabelColor(0);

    // Check that the file is readable
    if(!fin.good())
      {
        return;
      }

    // Read each line of the file separately
    for(unsigned int iLine=0;!fin.eof();iLine++)
      {
      // Read the line into a string
      std::getline(fin,line);

      // Check if the line is a comment or a blank line
      if(line[0] == '#' || line.length() == 0)
        continue;

      // Create a stream to parse that string
      std::istringstream iss(line);

      try 
        {
        // Read in the elements of the file
        uint16_t idx;
        int red, green, blue, visible, mesh;
        float alpha;
        iss >> idx;
        iss >> red;
        iss >> green;
        iss >> blue;
        iss >> alpha;
        iss >> visible;
        iss >> mesh;

        // Skip to a quotation mark
        iss.ignore(line.length(),'\"');

        // Allocate a label of appropriate size
        char *label = new char[line.length()+1];

        // Read the label
        iss.get(label,line.length(),'\"');

        if (idx > 0)
          {
          //std::cout << "label(" << idx << ") color(" 
          //   << +red << ',' << +green << ',' << +blue << ")" << std::endl;
          LabelColor lc(red, green, blue, alpha);
          m_LabelColorMap[idx] = lc;
          }

        /*
        // Create a new color label
        ColorLabel cl;

        // Store the results
        // cl.SetValid(true);
        cl.SetRGB(0,(unsigned char) red);
        cl.SetRGB(1,(unsigned char) green);
        cl.SetRGB(2,(unsigned char) blue);
        cl.SetAlpha( (unsigned char) (255 * alpha) );
        cl.SetVisible(visible != 0);
        cl.SetVisibleIn3D(mesh != 0);
        cl.SetLabel(label);
        */

        // Clean up the label
        delete[] label;
        }
      catch( std::exception e)
        {
        // Close the input stream
        fin.close();
        
        // create an exeption string
        std::cerr << "Exception when reading file: " << e.what() << std::endl;
        }
      }

    fin.close();
  }

  void GetDefaultLabelColor(uint16_t label, uint8_t &r, uint8_t &g, uint8_t &b)
  {
    parse_color(ColorList[(label-1) % ColorListSize], r, g, b);
    //std::cout << "Label[" << label << "] color=(" << +r << ',' << +g << ',' << +b << ")" << std::endl;
  }

  void GetLabelColor(uint16_t label, uint8_t &r, uint8_t &g, uint8_t &b, float &a)
  {
    if (m_LabelColorMap.count(label))
      {
      LabelColor lc = m_LabelColorMap[label];
      r = lc.r;
      g = lc.g;
      b = lc.b;
      a = lc.a;
      }
    else
      {
      GetDefaultLabelColor(label, r, g, b);
      a = 1.0;
      }
  }

private:
  LabelColorMap m_LabelColorMap;

  static void parse_color(const char* p, unsigned char& r, unsigned char& g, unsigned char& b)
  {
    assert(strlen(p) == 7); // Must be a seven-character string

    int val[6];
    for(int i = 0; i < 6; i++)
      {
      char c = p[i+1];
      if(c >= 'A' && c <= 'F')
        val[i] = 10 + (c - 'A');
      else if(c >= 'a' && c <= 'f')
        val[i] = 10 + (c - 'a');
      else if(c >= '0' && c <= '9')
        val[i] = c - '0';
      }

    r = 16 * val[0] + val[1];
    g = 16 * val[2] + val[3];
    b = 16 * val[4] + val[5];
  }
};