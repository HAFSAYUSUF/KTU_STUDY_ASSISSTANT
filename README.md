# 📘 KTU Study Assistant

KTU Study Assistant is an **AI-powered academic support application** designed to help **KTU students** learn more efficiently and manage their academic time effectively. The application combines multiple study tools into a single platform using **Python, Streamlit, and the Google Gemini 1.5 API**.

---

## 🚀 Features

- 🤖 **AI Chatbot Tutor**
  - Answers academic doubts in real time
  - Provides exam-oriented explanations, formulas, and summaries

- 🧠 **Quiz Generator**
  - Generates multiple-choice questions (MCQs)
  - Helps students evaluate their understanding

- 📄 **PDF Summarizer**
  - Summarizes uploaded lecture notes and study materials
  - Handles large PDFs using chunk-based processing

- 📑 **PDF Viewer**
  - View uploaded PDFs directly within the app

- 🎓 **SGPA & CGPA Calculator**
  - Calculates semester-wise SGPA
  - Computes cumulative CGPA across semesters

- 📅 **Exam Countdown & Study Planner**
  - Displays remaining days until exam
  - Helps divide syllabus across available days

- 📊 **Study Tracker**
  - Tracks daily study activity

---

## 🛠️ Tech Stack

- **Programming Language:** Python  
- **Frontend Framework:** Streamlit  
- **AI Model:** Google Gemini 1.5 API  
- **PDF Processing:** PyMuPDF  
- **Other Libraries:** `dotenv`, `datetime`, `re` (regex)

---

## 🧠 AI Integration

- Uses **Google Gemini 1.5** for:
  - Question answering
  - Quiz generation
  - PDF summarization
- **Prompt Engineering** ensures:
  - Syllabus-focused responses
  - Exam-oriented output
  - Clear and structured answers

---

## 🔐 Security

- API keys are stored using **environment variables**
- Sensitive credentials are not committed to GitHub
- Environment variables are accessed via Python’s `os` module

---

## ⚙️ How It Works

1. User selects a feature (Chatbot, Quiz, PDF Summarizer, etc.)
2. Input is structured using prompt engineering
3. Data is sent to the Gemini API
4. AI-generated output is rendered using Streamlit
5. Session data is maintained using `st.session_state`


