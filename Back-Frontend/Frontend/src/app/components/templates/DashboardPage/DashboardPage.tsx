import React from 'react';
import { DASHBOARD_ITEMS, SIDEBAR } from '../../../shared/constants';
import { IDashboardItem } from '../../../shared/model';
import Card from '../../atoms/Card';
import MainLayout from '../../layouts/MainLayout';
import Header from '../../organisms/Header';
import Sidebar from '../../organisms/Sidebar';
import './Dashboard.style.scss';
import arrowRight from '../../../../assets/images/arrow-right.svg';
import codeIcon from '../../../../assets/images/code-icon.svg';
import { DashboardSectionTitle } from './Dashboard.style';
import ModuleGradient from '../../atoms/ModuleGradient';
import docIcon from '../../../../assets/images/doc-icon.svg';
import Button from '../../atoms/Button';
import SearchBar from '../../atoms/SearchBar';
import { Link } from 'react-router-dom';

const BottomTip = ({ text }) => {
  return (
    <div className="bottom-text">
      {text}
      <img className="arrow-right" src={arrowRight} alt="" />
    </div>
  );
};

const DashboardItem = (props: {
  dashboardItem: IDashboardItem;
  index: number;
}) => {
  const { dashboardItem, index } = props;
  return (
    <Link to={dashboardItem.href} key={index}>
      <Card className="dashboard-module">
        <div className="dashboard-module__left">
          <p className="dashboard-module__title">{dashboardItem.title}</p>
          <p className="dashboard-module__subtitle">
            {dashboardItem.description}
          </p>
          <BottomTip text="Read More" />
        </div>
        <ModuleGradient className="dashboard-module__right">
          {dashboardItem.title}
        </ModuleGradient>
      </Card>
    </Link>
  );
};

const DashboardModules = () => {
  return (
    <React.Fragment>
      <DashboardSectionTitle>
        Modules
        <SearchBar />
      </DashboardSectionTitle>
      <div className="dashboard-modules">
        {DASHBOARD_ITEMS.map((dashboardItem, index) => (
          <DashboardItem dashboardItem={dashboardItem} index={index} />
        ))}
      </div>
    </React.Fragment>
  );
};

const ProjectSourceCode = () => {
  return (
    <React.Fragment>
      <DashboardSectionTitle>Project Source Code</DashboardSectionTitle>
      <Card className="project-source-code">
        <ModuleGradient className="project-source-code__left">
          <img alt="" src={codeIcon} />
        </ModuleGradient>
        <div className="project-source-code__right">
          <p className="dashboard-module__title">Source Code</p>
          <p className="dashboard-module__subtitle">
            Amet minim mollit non deserunt ullamco est sit aliqua dolor do amet
            sint. Velit officia consequat duis enim velit mollit. Exercitation
            veniam consequat sunt nostrud amet.
          </p>
          <BottomTip text="View Source Code" />
        </div>
      </Card>
    </React.Fragment>
  );
};

const Dashboard = (props) => {
  return (
    <React.Fragment>
      <Header pageTitle="Dashboard" />
      <div className="dashboard">
        <DashboardModules />
        <ProjectSourceCode />
      </div>
    </React.Fragment>
  );
};

const Docs = (props) => (
  <div className="dashboard-docs">
    <ModuleGradient>
      <div className="dashboard-docs__inner">
        <img alt="" src={docIcon} />
        <br />
        <br />
        <div className="dashboard-docs__text">
          <strong>Need help?</strong>
          <br />
          Please Check our docs
        </div>
        <br />
        <Button>Documentation</Button>
      </div>
    </ModuleGradient>
  </div>
);

const DashboardSidebar = (props) => (
  <Sidebar sidebarItems={SIDEBAR} activeItem={SIDEBAR[0]}>
    <Docs />
  </Sidebar>
);

class DashboardPage extends React.Component {
  render() {
    return (
      <MainLayout sidebar={<DashboardSidebar />} content={<Dashboard />} />
    );
  }
}

export default DashboardPage;
