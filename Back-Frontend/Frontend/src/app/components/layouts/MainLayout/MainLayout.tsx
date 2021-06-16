import React from "react";
import { MainLayoutContainer } from "./MainLayout.style";

const MainLayout = (props) => {
  const { sidebar, content } = props;
  return (
    <MainLayoutContainer>
      <section className="sidebar">{sidebar}</section>
      <section className="main-area">{content}</section>
    </MainLayoutContainer>
  );
};

export default MainLayout;
