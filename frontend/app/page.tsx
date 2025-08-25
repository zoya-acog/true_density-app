import Header from "@/components/header"; // Adjust path if needed
import InputBar from "@/components/input_bar"; // Adjusted import path to match typical Next.js structure


export default function Home() {
  return (
    <div>
      <Header />
      <main className="pt-15">
        <InputBar />
        
      </main>

    </div>
  );
}